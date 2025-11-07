library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ape)
library(data.table)
library(valr)
library(stringr)
library(data.table)
library(cowplot)
library(ggtree)


# This function intersects g2g long-read alignments with genome bins to extract bin-level coverage and identity 
# a miriad of bin metrics are calculated (see mutate cluster)
alignmentPartition <- function(df) {
  a2a_bed <- df %>% dplyr::select(REF,S1,E1,IDY,STRAIN) %>%
    dplyr::rename(chrom=REF,start=S1,end=E1,Identity=IDY)
  
  df_LRbins <- valr::bed_intersect(bins_1kb_CB_stripped,a2a_bed) %>%
    dplyr::rename(CHROM=chrom, START_BIN=start.x, END_BIN=end.x, 
                  start=start.y, stop=end.y, coverage=.overlap, Identity=Identity.y) %>%
    dplyr::group_by(CHROM,START_BIN,END_BIN)  %>% 
    dplyr::filter(coverage==max(coverage)) %>% 
    dplyr::filter(Identity==max(Identity)) %>% 
    dplyr::mutate(group_id=cur_group_id()) %>%
    dplyr::distinct(group_id,.keep_all = T) %>%
    dplyr::ungroup()
  
  df_LRbins_unclass <- bins_1kb_CB_stripped %>% 
    dplyr::rename(CHROM=chrom, START_BIN=start, END_BIN=end) %>%
    dplyr::left_join(df_LRbins,by=c("CHROM"="CHROM","START_BIN"="START_BIN","END_BIN"="END_BIN")) %>%
    dplyr::arrange(CHROM,START_BIN) %>% 
    dplyr::mutate(coverage=ifelse(is.na(coverage),0,coverage)) %>%
    dplyr::mutate(STRAIN.y=ifelse(is.na(STRAIN.y),unique(a2a_bed$STRAIN),STRAIN.y)) %>%
    dplyr::mutate(full_cov = ifelse(coverage == 1e3, T, F)) %>% 
    dplyr::rename(bin_IDY=Identity,bin_COV=coverage,STRAIN=STRAIN.y) %>%
    dplyr::select(CHROM,START_BIN,END_BIN,bin_IDY,bin_COV,STRAIN) 
  
  return(df_LRbins_unclass)
}

# given a set of coverage and identity thresholds, this function classifies bins as divergent or non-divergent
# divergent bins have sub-classifications based on the variant determinants 
classifyPartition <- function(df,cf,idy) {
  df_bins_clasif_LR <- df %>%
    dplyr::mutate(div=ifelse(bin_COV<(cf*10),"C",
                             ifelse(!(is.na(bin_IDY)) & bin_IDY<=idy,"I","nondiv"))) %>%
    dplyr::arrange(CHROM,START_BIN) %>%
    dplyr::group_by(CHROM) %>%
    dplyr::mutate(div_class = ifelse((lag(div) == "C" | lag(div) == "I") & (lead(div) == "C" | lead(div) == "I") & div=="nondiv","G", div)) %>%
    dplyr::mutate(div_class = ifelse(div_class=="C" & bin_COV == 0,"Z",div_class)) %>%
    dplyr::mutate(div_gf=ifelse(div_class=="nondiv","nondiv","div")) %>%
    dplyr::ungroup()
  
  return(df_bins_clasif_LR)
}

## this function estimates the frequency at which each bin is classified as hyperdivergent in the sampled population
## low frequency () bins are discarded to avoid super-clustering of divergent regions
frequencyEstimates <- function(df) {
  denom<-length(unique(df$STRAIN))
  df2 <- df %>% 
    dplyr::group_by(CHROM,START_BIN) %>%
    dplyr::mutate(binCt=stringr::str_count(div_gf,"nondiv")) %>%
    dplyr::mutate(binCt=ifelse(is.na(binCt),1,binCt)) %>%
    dplyr::mutate(binCt=ifelse(binCt==1,0,1)) %>%
    dplyr::mutate(binCtSum=sum(binCt)) %>%
    dplyr::mutate(binFreq=(binCtSum/denom)*100)
  return(df2)
}

## This function clusters contiguous div bins into div regions, 
## and assembles a divergent 'footprint' from the string of bin classifications
clusterBins <- function(df,mode) {
  
  if(mode=="SRF") {
    temp <- df %>% 
      dplyr::arrange(CHROM,START_BIN) %>%
      dplyr::mutate(div_gf=ifelse(div_gf =="div" & binFreq < 5,"nondiv",div_gf)) %>%
      dplyr::select(-binCt,-binCtSum,-binFreq)
  } else {
    temp <- df %>% 
      dplyr::arrange(CHROM,START_BIN) 
  }
  
  temp$enum <- sequence(rle(as.character(temp$div_gf))$lengths)
  
  div_bins <- temp %>%
    dplyr::arrange(CHROM,START_BIN) %>%
    dplyr::group_by(CHROM,data.table::rleid(div_gf)) %>%
    dplyr::mutate(gid=cur_group_id()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(CHROM) %>%
    dplyr::mutate(clusTips=ifelse(div_gf == "div" & lead(enum) == 2 & lead(div_gf)=="div","clust_start", 
                                  ifelse(div_gf == "div" & lead(enum) > 1 & lead(div_gf)=="div","clust_center",
                                         ifelse(div_gf == "div" & enum > 1 & lead(div_gf)=="nondiv","clust_end","unclust")))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!(clusTips=="unclust"))
  
  if (mode == "LR") {
    div_regions <- div_bins %>%
      dplyr::group_by(gid) %>%
      dplyr::mutate(prop_covz=sum(stringr::str_detect(div_class,"Z"))/n()) %>%
      dplyr::mutate(prop_lowcov=sum(stringr::str_detect(div_class,"C"))/n()) %>%
      dplyr::mutate(prop_gf=sum(stringr::str_detect(div_class,"G"))/n()) %>%
      dplyr::mutate(prop_idy=sum(stringr::str_detect(div_class,"I"))/n()) %>%
      dplyr::mutate(
        meanIDY = if (all(is.na(bin_IDY))) NA_real_ else
          sum(bin_IDY * bin_COV, na.rm = TRUE) / sum(bin_COV[!is.na(bin_IDY)])) %>%
      dplyr::mutate(meanCov=mean(bin_COV),minStart=min(START_BIN),maxEnd=max(END_BIN),divSize=maxEnd-minStart,group_size=n()) %>%
      dplyr::mutate(bin_foot=paste0(div_class,collapse = "")) %>%
      dplyr::distinct(gid,.keep_all = T) %>%
      dplyr::ungroup() %>%
      dplyr::select(CHROM,minStart,maxEnd,divSize,meanIDY,meanCov,prop_covz,prop_lowcov,prop_gf,prop_idy,bin_foot,group_size)
  } else {
    div_regions <- div_bins %>%
      dplyr::group_by(gid) %>%
      dplyr::mutate(prop_covz=sum(stringr::str_detect(div_class,"Z"))/n()) %>%
      dplyr::mutate(prop_lowcov=sum(stringr::str_detect(div_class,"C"))/n()) %>%
      dplyr::mutate(prop_gf=sum(stringr::str_detect(div_class,"G"))/n()) %>%
      dplyr::mutate(prop_idy=sum(stringr::str_detect(div_class,"I"))/n()) %>%
      dplyr::mutate(meanVC=mean(COUNT),meanCF=mean(pc1X),minStart=min(START_BIN),maxEnd=max(END_BIN),divSize=maxEnd-minStart,group_size=n()) %>%
      dplyr::mutate(bin_foot=paste0(div_class,collapse = "")) %>%
      dplyr::distinct(gid,.keep_all = T) %>%
      dplyr::ungroup() %>%
      dplyr::select(CHROM,minStart,maxEnd,divSize,meanVC,meanCF,prop_covz,prop_lowcov,prop_gf,prop_idy,bin_foot,group_size)
  }
  
  out <- list(div_bins,div_regions)
  return(out)
}

#get nucmer g2g coords
transformed_coords <- readr::read_tsv("../processed_data/Tropical_hifi_coords.tsv",
                                      col_names = c("S1","E1","S2","E2","L1","L2","IDY","LENR","LENQ","REF","HIFI","STRAIN")) %>%
  dplyr::filter(STRAIN!= "ECA3733" & STRAIN!="BRC20341") # no SR data or wrong RG

#list of strains with long-read genomes
div_str <- unique(transformed_coords$STRAIN)

#get 1kb windows (bedtools makewindows)
bins_1kb_CB_stripped <- readr::read_tsv("../processed_data/QX1410_genomic_windows.1kb.bed",col_names = F) %>%
  dplyr::rename(chrom=X1,start=X2,end=X3) 

#Filter out alignments smaller than 1kb
#tcords is now a list of datframes that contain the g2g alignment coordinates for each strain
tcords <- transformed_coords %>%
  dplyr::filter(L2>1000) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::group_split()

#use alignementPartition() to extract coverage and idy stats from g2g coordinate sets for each bin
strain_df <- list()
for (i in 1:length(tcords)) {
  print(i)
  tmp <- alignmentPartition(tcords[[i]])
  strain_df[[i]]  <- tmp
}
all_LR_stats <-ldply(strain_df,data.frame)

####### CALL LR BASED HDRs ##########
#set a range of identities to test HDR calling using LR data
cf <- 60
idthresh <- c(seq(90,99,1))

#iterate through identity range and classify partitions previously generated
strain_dv <- list()
clasiBins <- list()
for (i in 1:length(strain_df)) {
  print(i)
  strain_varID <- list()
  for (k in 1:length(idthresh)) {
    strain_varID[[k]] <- classifyPartition(strain_df[[i]],cf,idthresh[k])
    strain_varID[[k]]$threshIDY <- idthresh[k]
  }
  strain_dv[[i]] <- strain_varID
  clasiBins[[i]] <- ldply(strain_varID,data.frame)
}

# this df holds bin data after HDR classification at 95% IDY, mainly for diagnostic purposes, can be omitted
#all_bins_LR <- ldply(clasiBins,data.frame) 

#iterate through bin classifications at various idy thresholds to cluster bins into regions using clusterBins()
dvReg <- list()
dvReg_part <- list()
for (i in 1:length(strain_df)) {
  print(i)
  strain_declust <- data.frame(CHROM=character(),minStart=integer(),maxEnd=integer(),divSize=integer(),meanIDY=numeric(),treshIDY=integer())
  bin_declust <- data.frame(CHROM=character(),START_BIN=integer(),END_BIN=integer(),bin_IDY=double(),bin_COV=double(),group_size=integer(),div_class=character(),gid=integer(),clusTips=character())
  strain <- unique(strain_dv[[i]][[1]]$STRAIN)
  
  for (k in 1:length(idthresh)) {
    temp <- clusterBins(strain_dv[[i]][[k]],"LR")
    temp_bins <- temp[[1]] %>%
      dplyr::select(CHROM,START_BIN,END_BIN,bin_IDY,bin_COV,div_class,gid,clusTips)
    temp_div <- temp[[2]]
    temp_div$threshIDY <- idthresh[k]
    temp_bins$threshIDY <- idthresh[k]
    strain_declust <- rbind(strain_declust,temp_div)
    bin_declust <- rbind(bin_declust,temp_bins)
  }
  strain_declust$STRAIN <- strain
  bin_declust$STRAIN <- strain
  
  dvReg[[i]] <- strain_declust %>%
    dplyr::mutate(regionClass=ifelse(prop_covz>0.95,"coverage gap",
                                     ifelse(prop_idy>0.5, "identity call", "coverage call")))
  dvReg_part[[i]] <- bin_declust
}
all_calls_LR <- ldply(dvReg, data.frame) 


LR_summ <- ggplot(all_calls_LR %>% dplyr::filter((threshIDY==94|threshIDY==95|threshIDY==96|threshIDY==97) & CHROM=="II" & (grepl("QG",STRAIN)))) + 
  geom_rect(aes(xmin=minStart/1e6,xmax=maxEnd/1e6,ymin=threshIDY-0.4,ymax=threshIDY+0.4,fill=regionClass)) + facet_grid(STRAIN~CHROM, scales = 'free') +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        strip.text = element_text(size=10)) +
  scale_x_continuous(expand=c(0,0)) +
  coord_cartesian(xlim=c(12,14))
LR_summ


s1 <- ggplot(all_calls_LR %>% dplyr::filter(divSize > 200e3)) + 
  geom_histogram(aes(x=divSize/1e3),binwidth = 5) + 
  facet_grid(CHROM~threshIDY, scales = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45,hjust = 1)) +
  xlab("HDR size (kb)") +
  scale_y_continuous(breaks = seq(0, max(all_calls_LR$divSize/1e3, na.rm = TRUE), by = 3))

s2 <- ggplot(all_calls_LR %>% dplyr::filter(divSize > 50e3 & divSize <=200e3)) + 
  geom_histogram(aes(x=divSize/1e3),binwidth = 5) + 
  facet_grid(CHROM~threshIDY, scales = "free_x")  +
  theme_classic() +
  theme(strip.text.x=element_blank(),
        axis.text.x = element_text(angle=45,hjust = 1)) +
  xlab("HDR size (kb)")

s3 <- ggplot(all_calls_LR %>% dplyr::filter(divSize >= 5e3 & divSize <=50e3)) + 
  geom_histogram(aes(x=divSize/1e3),binwidth = 5) + 
  facet_grid(CHROM~threshIDY, scales = "free_x") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45,hjust = 1))


idy_bins <- all_LR_stats %>% dplyr::filter(!is.na(bin_IDY) & bin_COV >= 600) %>% dplyr::filter(STRAIN!= "ECA3733" & STRAIN!="BRC20341")

strain_thresh <- idy_bins %>%
  dplyr::filter(bin_COV>=600) %>%
  dplyr::group_by(STRAIN, CHROM) %>%
  dplyr::summarise(
    sample_mean = mean(bin_IDY, na.rm = TRUE),
    sample_sd   = sd(bin_IDY, na.rm = TRUE),
    thresh      = sample_mean - 2 * sample_sd,
    upper_q95   = quantile(bin_IDY, probs = 0.05, na.rm = TRUE),
    .groups = "drop")

strain_thresh_all <- idy_bins %>%
  dplyr::filter(bin_COV>=600) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::summarise(
    sample_mean = mean(bin_IDY, na.rm = TRUE),
    sample_sd   = sd(bin_IDY, na.rm = TRUE),
    thresh      = sample_mean - 2 * sample_sd,
    upper_q95   = quantile(bin_IDY, probs = 0.05, na.rm = TRUE),
    .groups = "drop") %>%
  dplyr::mutate(CHROM="All")

s4 <- ggplot(
  rbind(strain_thresh, strain_thresh_all),
  aes(x = STRAIN, y = thresh, group = CHROM, color = CHROM)
) +
  geom_point() +
  geom_line() +
  theme_classic() +
  theme(axis.text = element_text(angle = 45, hjust = 1)) +
  ylab("Bin identity\nlower 5% quantile") +
  xlab("Strain") +
  scale_color_manual(
    values = c("I"   = "#1b9e77",
               "II"  = "#d95f02",
               "III" = "#7570b3",
               "IV"  = "#e7298a",
               "V"   = "#66a61e",
               "X"   = "#e6ab02",
               "All" = "black")) +
  geom_hline(yintercept = 96,linetype="dashed")

pc1<- cowplot::plot_grid(s1,s4,ncol=2, rel_widths = c(0.5,1),align = "v",axis = "l",labels = c("c","d"))
pc2<- cowplot::plot_grid(s3,s2,ncol=1, rel_heights = c(1.1,1),align = "v",axis = "lr",labels=c("a","b"))

comp_reg<-cowplot::plot_grid(pc2,pc1,ncol=1,nrow=2,rel_heights = c(2,1),align = "hv",axis = "rlbt")

ggsave(plot = comp_reg, filename = "../figures/FigureS31_LR_STATS_95idy_20251023.png",width = 7.5,height = 9,dpi = 600,device = 'png')


####### CALL SR BASED HDRS (FOR LR STRAINS) ##########
covfrac <- c(seq(0.05,0.9,0.05))
varct <- c(seq(5,30,1))

#all strain coverage data
coverage_df <- readr::read_table("../processed_data/Tropical.thresh_cov.tsv",col_names = c("CHROM","START_BIN","END_BIN","NAME","c1X","c2X","c5X","STRAIN"))
#all strain variant count data
varct_df <- readr::read_table("../processed_data/Tropical.variant_counts.tsv", col_names = c("CHROM","START_BIN","END_BIN","COUNT","STRAIN"))

#join coverage and variant count data
SR_stats_all <- varct_df %>%
  dplyr::left_join(coverage_df,by=c("STRAIN","CHROM","START_BIN","END_BIN")) %>%
  dplyr::select(-NAME) %>%
  dplyr::filter(!(CHROM=="MtDNA")) %>%
  dplyr::mutate(c1X=ifelse(is.na(c1X),0,c1X)) %>%
  dplyr::mutate(c2X=ifelse(is.na(c2X),0,c2X)) %>%
  dplyr::mutate(c5X=ifelse(is.na(c5X),0,c5X)) %>%
  dplyr::mutate(pc1X=c1X/1e3,pc2X=c2X/1e3,pc5X=c5X/1e3) %>%
  dplyr::filter(STRAIN %in% div_str) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(STRAIN_count=sum(COUNT)) %>% 
  dplyr::ungroup() 
  

#iterate through range of coverage and variant count thresholds
#at each iteration, classify bins and cluster adjacent bins
SR_list <- list()
ct=1
for (i in 1:length(covfrac)) {
  for (k in 1:length(varct)){
    print(paste0("cf:",covfrac[i]," / vc:",varct[k]))
    #classify
    all_stats <- SR_stats_all %>%
      dplyr::mutate(div_class=ifelse(COUNT >= varct[k],"I", 
                                     ifelse(pc1X < covfrac[i],"C","R"))) %>% 
      dplyr::mutate(div=ifelse(div_class=="C" | div_class=="I","div","nondiv")) %>%               
      dplyr::group_by(STRAIN,CHROM) %>%
      dplyr::mutate(div_gf=ifelse(div=="nondiv" & lead(div)=="div" & lag(div)=="div","div",div)) %>%
      dplyr::mutate(div_class=ifelse(div==div_gf,div_class,"G")) %>%
      dplyr::ungroup() 
    
    #cluster
    regList <- list()
    #binList <- list()
    for (j in 1:length(div_str)) {
      temp <- all_stats %>% dplyr::filter(STRAIN==div_str[j])
      div_call <- clusterBins(temp,"SR")
      # div_call[[1]]$STRAIN <- div_str[[j]]
      # div_call[[1]]$CFT <- covfrac[[i]]
      # div_call[[1]]$VCT <- varct[[k]]
      div_call[[2]]$STRAIN <- div_str[[j]]
      div_call[[2]]$CFT <- covfrac[[i]]
      div_call[[2]]$VCT <- varct[[k]]
      #binList[[i]] <- div_call[[1]]
      regList[[j]] <- div_call[[2]]
    }
    
    SR_list[[ct]] <- ldply(regList, data.frame) 
    ct = ct + 1
  }
}
all_SR_meta_calls <- ldply(SR_list,data.frame)

#save.image(file="/vast/eande106/projects/Nicolas/hyperdivergent_regions/briggsae/multi_reference/QX1410/HDR_20250730_analysis_CB_chkpt1_TROP.Rda")
#based on Daehan's paper and results from "pcomp" plot above we use 95% IDY LR calls as the truth set
t95_LR_calls <- all_calls_LR %>%
  dplyr::filter(threshIDY==96) %>%
  dplyr::select(CHROM,minStart,maxEnd,STRAIN,regionClass) %>%
  dplyr::rename(chrom=CHROM,start=minStart,end=maxEnd)

#iterate through range of coverage and variant count thresholds
#at each iteration, we identify intersects between LR and SR HDR calls
#we estimate overlap fraction, excess fraction, recall, and precision
#we correct for 1:Many overlaps by grouping intersects by the LR call, and aggregating their overlap and excess fractions
#for Many:1, we keep the SR call with the longest overlap
interList_JC <- list()
ct=1
for (i in 1:length(covfrac)) {
  for (k in 1:length(varct)){
    print(paste0("cf:",covfrac[i]," / vc:",varct[k]))
    temp <- all_SR_meta_calls %>%
      dplyr::filter(VCT==varct[k] & CFT==covfrac[i])
    
    for (j in 1:length(div_str)) {
      #get SR calls
      temp_SR <- temp %>%
        dplyr::filter(STRAIN==div_str[j]) %>%
        dplyr::select(CHROM, minStart,maxEnd,STRAIN,VCT,CFT) %>%
        dplyr::rename(chrom=CHROM,start=minStart,end=maxEnd) %>%
        dplyr::filter(end-start >= 5e3)
      
      #get LR calls
      temp_LR <- t95_LR_calls  %>% 
        dplyr::filter(STRAIN==div_str[j]) %>%
        dplyr::filter(end-start >= 5e3)
      
      #call intersects
      interSet <- valr::bed_intersect(temp_SR,temp_LR) %>%
        dplyr::mutate(LR_size=end.y-start.y) %>%
        dplyr::mutate(SR_size=end.x-start.x) %>%
        #dplyr::mutate(LR_centroid=start.x+(LR_size/2)) %>%
        #dplyr::mutate(SR_centroid=start.y+(SR_size/2)) %>%
        dplyr::mutate(overlap_fraction=.overlap/LR_size) %>% #overlap fraction (OF)
        #dplyr::mutate(centroid_shift=abs(LR_centroid-SR_centroid)) %>%
        #dplyr::mutate(sizeDiff=abs(LR_size-SR_size)) %>%
        dplyr::mutate(exc_fraction=(SR_size-.overlap)/SR_size) %>% # excess fraction (EF)
        dplyr::mutate(LRID=paste(chrom,start.y,end.y,sep = "_")) %>% #generate a unique ID for each LR call
        dplyr::mutate(SRID=paste(chrom,start.x,end.x,sep = "_")) %>% #generate a unique ID for each SR call
        dplyr::group_by(LRID) %>% #group by LR call
        dplyr::mutate(LRgsize=n(),
                      agg_OF=sum(.overlap)/LR_size,
                      agg_EF=(sum(SR_size)-sum(.overlap))/sum(SR_size),
                      locmax_OF=ifelse(.overlap==max(.overlap),T,F)) %>% #correct OF and EF, flag longest overlap
        dplyr::ungroup() %>%
        dplyr::filter(locmax_OF==T) %>% #keep one entry per group (longest overlap)
        dplyr::group_by(SRID) %>% #group by SR call
        dplyr::mutate(SRgsize=n(), 
                      locmax_SRo=ifelse(.overlap==max(.overlap),T,F)) %>% #flag longest overlap
        dplyr::ungroup() %>%
        dplyr::filter(locmax_SRo==T) #keep one entry per group (longest overlap)
      
      interSet$noverRec <- nrow(interSet)/nrow(temp_LR) #recall
      interSet$noverPre <- nrow(interSet)/nrow(temp_SR) #precision
      interSet$noverF1 <- 2*((nrow(interSet)/nrow(temp_SR) * nrow(interSet)/nrow(temp_LR))/(nrow(interSet)/nrow(temp_SR) + nrow(interSet)/nrow(temp_LR)))
        
      interList_JC[[ct]] <-interSet
      
      ct=ct+1
    }
  }
}

#aggregate all intersect data
all_intersects <- ldply(interList_JC,data.frame) %>%
  dplyr::mutate(PP=paste0(VCT.x,":",CFT.x)) #%>% #set parameter pair (PP) ID
  #dplyr::filter(overlap_fraction > 0) #placeholder, could be used to filter low quality overlaps

############ PARAMETER PAIR SELECTION ##########

statSumm <- all_intersects %>%
  dplyr::group_by(STRAIN.x,PP) %>%
  dplyr::mutate(meanOF=mean(agg_OF)) %>%
  dplyr::mutate(meanEF=mean(agg_EF)) %>%
  dplyr::mutate(gid=cur_group_id()) %>%
  dplyr::ungroup()

maxOF <- statSumm %>%
  dplyr::group_by(gid) %>%
  dplyr::arrange(desc(meanOF)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(gid,.keep_all = T)

MOF <- ggplot(maxOF) + geom_point(aes(x=VCT.x,y=meanOF,color=as.factor(CFT.x*100)),size=0.5)  +
  facet_wrap(~STRAIN.x) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill=NA)) +
  guides(color=guide_legend(title="Percent\nbases\ncovered")) +
  xlab("Variant count") +
  ylab("Mean overlap fraction") +
  geom_line(aes(x=VCT.x,y=meanOF,color=as.factor(CFT.x*100),group=CFT.x))

#meanEF
MEF <- ggplot(maxOF) + geom_point(aes(x=VCT.x,y=meanEF,color=as.factor(CFT.x*100)),size=0.5)  +
  facet_wrap(~STRAIN.x) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill=NA)) +
  guides(color=guide_legend(title="Percent\nbases\ncovered")) +
  xlab("Variant count") +
  ylab("Mean excess fraction") +
  geom_line(aes(x=VCT.x,y=meanEF,color=as.factor(CFT.x*100),group=CFT.x)) 

REC <- ggplot(maxOF) + geom_point(aes(x=VCT.x,y=noverRec,color=as.factor(CFT.x*100)),size=0.5)  +
  facet_wrap(~STRAIN.x) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill=NA)) +
  guides(color=guide_legend(title="Percent\nbases\ncovered")) +
  xlab("Variant count") +
  ylab("Recall") +
  geom_line(aes(x=VCT.x,y=noverRec,color=as.factor(CFT.x*100),group=CFT.x)) 

PRE <- ggplot(maxOF) + geom_point(aes(x=VCT.x,y=noverPre,color=as.factor(CFT.x*100)),size=0.5)  +
  facet_wrap(~STRAIN.x) +
  #geom_vline(xintercept = 14) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill=NA)) +
  guides(color=guide_legend(title="Percent\nbases\ncovered")) +
  xlab("Variant count") +
  ylab("Precision") +
  geom_line(aes(x=VCT.x,y=noverPre,color=as.factor(CFT.x*100),group=CFT.x)) 

F1 <- ggplot(maxOF) + geom_point(aes(x=VCT.x,y=noverF1,color=as.factor(CFT.x*100)),size=0.5)  +
  facet_wrap(~STRAIN.x) +
  #geom_vline(xintercept = 14) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill=NA)) +
  guides(color=guide_legend(title="Percent\nbases\ncovered")) +
  xlab("Variant count") +
  ylab("F1") +
  geom_line(aes(x=VCT.x,y=noverF1,color=as.factor(CFT.x*100),group=CFT.x)) 

topHits <- statSumm %>%
  dplyr::group_by(gid) %>%
  dplyr::arrange(desc(noverF1)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(gid,.keep_all = T) %>%
  dplyr::group_by(STRAIN.x) %>%
  dplyr::arrange(desc(noverF1)) %>%
  dplyr::slice_max(order_by = noverF1, n = 10) %>% #n is tweaked manually to find consensus optimal

bestHit <- topHits %>%
  dplyr::group_by(PP) %>%
  dplyr::mutate(freqPP=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(maxfreq=max(freqPP)) %>%
  dplyr::filter(freqPP==max(freqPP)) %>%
  dplyr::mutate(class="HIFREQ") %>%
  dplyr::select(-maxfreq)

worseHits <- topHits %>%
  dplyr::filter(!(gid %in% bestHit$gid)) %>%
  dplyr::mutate(class="TOPHIT") %>%
  dplyr::group_by(STRAIN.x) %>%
  dplyr::mutate(class=ifelse(noverF1==max(noverF1),"Strain\noptimal",class)) %>%
  dplyr::ungroup()

BEST <- ggplot() + 
  geom_point(data = worseHits %>% filter(class == "Strain\noptimal"),
             aes(x = VCT.x, y = CFT.x, color = meanOF, shape = class), size = 2) +
  geom_point(data = bestHit %>% mutate(class = "Consensus\noptimal\n"),
             aes(x = VCT.x, y = CFT.x, color = meanOF, shape = class), size = 2) +
  facet_wrap(~STRAIN.x) +
  xlab("Variant count") + ylab("Percent bases covered") +
  scale_color_gradient(name = "Mean\noverlap\nfraction", low = "orange", high = "blue") +
  scale_shape_manual(values = c("Strain\noptimal" = 15, "Consensus\noptimal\n" = 17)) +
  guides(shape = guide_legend(title = NULL)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA))

####################### CALL SR HDRS SPECIES-WIDE ###############################
#get list of WI (generated from VCF)
strainL <- readLines("../processed_data/Tropical_samples.txt")

COV_thresh = 0.9
VC_thresh = 11

#merge bin-based variant counts and coverage data
SR_stats_WI <- varct_df %>%
  dplyr::left_join(coverage_df,by=c("STRAIN","CHROM","START_BIN","END_BIN")) %>%
  dplyr::select(-NAME) %>%
  dplyr::filter(!(CHROM=="MtDNA")) %>%
  dplyr::mutate(c1X=ifelse(is.na(c1X),0,c1X)) %>%
  dplyr::mutate(c2X=ifelse(is.na(c2X),0,c2X)) %>%
  dplyr::mutate(c5X=ifelse(is.na(c5X),0,c5X)) %>%
  dplyr::mutate(pc1X=c1X/1e3,pc2X=c2X/1e3,pc5X=c5X/1e3) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(STRAIN_count=sum(COUNT)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(COUNT_ADJ = COUNT / STRAIN_count) 

#classify bins
all_stats <- SR_stats_WI %>%
  dplyr::mutate(div_class=ifelse(COUNT >= VC_thresh,"I", #over variant count thresh
                                 ifelse(pc1X < COV_thresh,"C","R"))) %>%  #under coverage thresh
  dplyr::mutate(div=ifelse(div_class=="C" | div_class=="I","div","nondiv")) %>%               
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::mutate(div_gf=ifelse(div=="nondiv" & lead(div)=="div" & lag(div)=="div","div",div)) %>% #fill 1kb gaps
  dplyr::mutate(div_class=ifelse(div==div_gf,div_class,"G")) %>% #assign a classification to gaps
  dplyr::ungroup()

#cluster bins (no frequency filter)
regList <- list()
binList <- list()
for (i in 1:length(strainL)) {
  print((i/length(strainL))*100)
  temp <- all_stats %>% dplyr::filter(STRAIN==strainL[i])
  div_call <- clusterBins(temp,"SR")
  #div_call[[1]]$STRAIN <- strainL[[i]]
  div_call[[2]]$STRAIN <- strainL[i]
  #binList[[i]] <- div_call[[1]]
  regList[[i]] <- div_call[[2]]
}

all_calls_SR <- ldply(regList, data.frame)

#QX calls are spurious (repeats) and should be removed from other strain calls
QX_calls <- all_calls_SR %>% dplyr::filter(STRAIN=="QX1410") %>% dplyr::select(CHROM,minStart,maxEnd)

#find overlaps and remove
dt_all <- as.data.table(all_calls_SR[all_calls_SR$STRAIN != "QX1410", ]) #all other strain calls
dt_qx  <- as.data.table(QX_calls) #qx calls
setnames(dt_qx, c("CHROM", "minStart", "maxEnd"), c("CHROM", "start", "end")) #rename keys
setnames(dt_all, c("minStart", "maxEnd"), c("start", "end")) #rename keys
setkey(dt_qx, CHROM, start, end) #set the keys for matching
setkey(dt_all, CHROM, start, end)  #set the keys for matching
overlaps <- foverlaps(dt_all, dt_qx, nomatch = 0L, type = "any") #find any overlap
to_drop <- unique(overlaps[, .(CHROM, start = i.start, end = i.end, STRAIN)]) #grab all other strain calls that overlap
key_cols <- c("CHROM", "start", "end", "STRAIN")
dt_filtered <- dt_all[!to_drop, on = key_cols] #remove them
setnames(dt_filtered, c("start", "end"), c("minStart", "maxEnd")) 
all_calls_SR_noQX <- as.data.frame(dt_filtered) #convert back to df for next steps

#flag adjacent regions that are 5kb apart
gap_clust <- all_calls_SR_noQX %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::mutate(forGapSize=lead(minStart)-maxEnd) %>%
  dplyr::mutate(flag3g=ifelse(forGapSize<=5000,"clust","noclust")) %>%
  dplyr::mutate(dec3g=ifelse(flag3g=="clust" ,"join",
                             ifelse(flag3g=="noclust" & lag(flag3g)=="clust","join","nojoin"))) %>%
  dplyr::mutate(dec3g=ifelse(is.na(dec3g),"nojoin",dec3g)) %>%
  dplyr::ungroup()

#get flagged and merged them
joinClust <- gap_clust %>% 
  dplyr::filter(dec3g=="join") %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::mutate(segbreak=ifelse(flag3g=="noclust",paste0(dec3g,row_number()),NA)) %>%
  tidyr::fill(segbreak,.direction = 'up') %>%
  dplyr::mutate(gid=data.table::rleid(segbreak)) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(conID=paste0(CHROM,"-",STRAIN,"-",gid)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(conID) %>%
  #dplyr::mutate(gid2=cur_group_id()) #%>%
  dplyr::mutate(newStart=min(minStart),newEnd=max(maxEnd)) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(gapFoot=paste0(rep("R",forGapSize/1e3),collapse = "")) %>%
  dplyr::mutate(new_foot=ifelse(flag3g=="clust",paste0(bin_foot,gapFoot),bin_foot)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(conID) %>%
  dplyr::mutate(clust_foot=paste0(new_foot,collapse = "")) %>%
  dplyr::mutate(newDivSize=newEnd-newStart) %>%
  dplyr::mutate(newMeanVC=mean(meanVC)) %>%
  dplyr::mutate(newMeanCF=mean(meanCF)) %>%
  dplyr::mutate(nclust=n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(conID,.keep_all = T) %>%
  dplyr::select(-minStart,-maxEnd,-divSize,-meanVC,-meanCF,-bin_foot) %>%
  dplyr::rename(minStart=newStart,maxEnd=newEnd,divSize=newDivSize,meanVC=newMeanVC,meanCF=newMeanCF,bin_foot=clust_foot) %>%
  dplyr::select(CHROM,minStart,maxEnd,divSize,meanVC,meanCF,bin_foot,STRAIN,nclust)

#keep unflagged 
nojoin <- gap_clust %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::filter(!(dec3g=="join")) %>%
  dplyr::ungroup() %>%
  dplyr::select(CHROM,minStart,maxEnd,divSize,meanVC,meanCF,bin_foot,STRAIN) %>%
  dplyr::mutate(nclust=1)

#bind unflagged and merged calls, arrange strains by number of HDRs
all_calls_SR_clustered <- rbind(joinClust,nojoin) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ncalls=n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(ncalls),STRAIN,CHROM,minStart) %>%
  dplyr::mutate(sorter=paste0(ncalls,STRAIN)) %>%
  dplyr::mutate(rleID=data.table::rleid(sorter)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ystrain=cur_group_id()) %>%
  dplyr::ungroup() 

all_calls_SR_clustered_sfilt <- all_calls_SR_clustered %>%
  dplyr::filter(divSize >= 5e3)

############# PLOT and WRITE
# 
# p1 <- ggplot(all_calls_SR_clustered %>% dplyr::mutate(REF="QX1410") %>% dplyr::filter(divSize >= 5e3 & !(CHROM=="MtDNA"))) + 
#   geom_rect(aes(xmin=minStart/1e6,xmax=maxEnd/1e6,ymin=rleID-0.45,ymax=rleID+0.45)) + 
#   facet_grid(REF~CHROM,scales = 'free') + 
#   theme(panel.background = element_blank(),
#         panel.grid = element_blank(),
#         axis.line.x = element_line(),
#         panel.border = element_rect(fill = NA),
#         axis.ticks.y = element_blank(),
#         axis.text.y=element_blank())  +
#   #xlab("Physical position (Mb)") +
#   ylab("508 Tropical isotype strains") +
#   #scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_x_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 5),expand = c(0, 0))

ggsave(plot = MOF, filename = "../figures/FigureS32_MOF_HDR_CB_20250730.png",width = 7.5,height = 6.5,dpi = 600,device = 'png')
ggsave(plot = MEF, filename = "../figures/FigureS33_MEF_HDR_CB_20250730.png",width = 7.5,height = 6.5,dpi = 600,device = 'png')
ggsave(plot = REC, filename = "../figures/FigureS34_REC_HDR_CB_20250730.png",width = 7.5,height = 6.5,dpi = 600,device = 'png')
ggsave(plot = PRE, filename = "../figures/FigureS35_PRE_HDR_CB_20250730.png",width = 7.5,height = 6.5,dpi = 600,device = 'png')
ggsave(plot = F1, filename = "../figures/FigureS36_F1_HDR_CB_20250730.png",width = 7.5,height = 6.5,dpi = 600,device = 'png')
ggsave(plot = BEST, filename = "../figures/FigureS37_BEST_HDR_CB_20250730.png",width = 7.5,height = 6.5,dpi = 600,device = 'png')

lineages <- readr::read_tsv("../processed_data/isotype_byLineage_GeoLocAdmCol_20250909.tsv") %>%
  dplyr::mutate(sublineage_color=ifelse(Sublineage=="TC","#ff0000",sublineage_color)) %>%
  dplyr::mutate(Sublineage=ifelse(Sublineage=="TC","TT",Sublineage))  %>% 
  dplyr::mutate(REF=ifelse(Lineage=="Tropical","QX1410",
                           ifelse(Lineage=="Temperate","JU2536",
                                  ifelse(Lineage=="TH","NIC1660",
                                         ifelse(Lineage=="AD","ECA2670",
                                                ifelse(Lineage=="KD","JU1348",
                                                       ifelse(Lineage=="TD1","BRC20530",
                                                              ifelse(Lineage=="TD2","BRC20492",NA))))))))

coverage_df_NR <- readr::read_table("../processed_data/Other_RG.thresh_cov.tsv",col_names = c("CHROM","START_BIN","END_BIN","NAME","c1X","c2X","c5X","STRAIN")) %>% dplyr::filter(!grepl("ptg",CHROM))
varct_df_NR <- readr::read_table("../processed_data/Other_RG.variant_counts.tsv", col_names = c("CHROM","START_BIN","END_BIN","COUNT","STRAIN")) %>% dplyr::filter(!grepl("ptg",CHROM))


coverage_df_NR2 <- coverage_df_NR %>%
  dplyr::left_join(lineages,by=c("STRAIN"="isotype")) %>%
  dplyr::filter(STRAIN!=REF)

varct_df_NR2 <- varct_df_NR %>%
  dplyr::left_join(lineages,by=c("STRAIN"="isotype")) %>%
  dplyr::filter(STRAIN!=REF)

SR_stats_WI_NR <- varct_df_NR2 %>%
  dplyr::left_join(coverage_df_NR2,by=c("REF","STRAIN","CHROM","START_BIN","END_BIN")) %>%
  dplyr::select(-NAME) %>%
  dplyr::filter(!(CHROM=="MtDNA")) %>%
  dplyr::mutate(c1X=ifelse(is.na(c1X),0,c1X)) %>%
  dplyr::mutate(c2X=ifelse(is.na(c2X),0,c2X)) %>%
  dplyr::mutate(c5X=ifelse(is.na(c5X),0,c5X)) %>%
  dplyr::mutate(pc1X=c1X/1e3,pc2X=c2X/1e3,pc5X=c5X/1e3) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(STRAIN_count=sum(COUNT)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(COUNT_ADJ = COUNT / STRAIN_count) 

plot_dist <- SR_stats_WI_NR %>% 
  #dplyr::mutate(COUNT=ifelse(c1X==0,NA,COUNT)) %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::summarise(p5 = quantile(COUNT, 0.05, na.rm = TRUE),
                   median = median(COUNT, na.rm = TRUE),
                   p95 = quantile(COUNT, 0.95, na.rm = TRUE),
                   p97 = quantile(COUNT, 0.97, na.rm = TRUE),
                   p99 = quantile(COUNT, 0.99, na.rm = TRUE),
                   REF=dplyr::first(REF)) %>% 
  dplyr::left_join(lineages %>% dplyr::select(isotype,Lineage),by=c("REF"="isotype")) %>%
  dplyr::mutate(Lineage=ifelse(REF=="NIC1667","TH (2)",Lineage)) %>%
  dplyr::ungroup()


plot_dist_trop <- SR_stats_WI %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::summarise(p5 = quantile(COUNT, 0.05, na.rm = TRUE),
                   median = median(COUNT, na.rm = TRUE),
                   p95 = quantile(COUNT, 0.95, na.rm = TRUE),
                   p97 = quantile(COUNT, 0.97, na.rm = TRUE),
                   p99 = quantile(COUNT, 0.99, na.rm = TRUE)) %>% 
  dplyr::left_join(lineages %>% dplyr::select(isotype,Lineage,Sublineage),by=c("STRAIN"="isotype"))  %>%
  dplyr::ungroup() %>%
  dplyr::filter(STRAIN!="QX1410")

otherlin<-ggplot(data=plot_dist,aes(x=STRAIN,y=p95,color=CHROM,group=CHROM)) + 
  geom_point(size=0.3) + 
  geom_line(linewidth=0.3) + 
  theme_classic() + 
  theme(legend.position = "right",
        axis.text.x = element_blank()) + 
  facet_wrap(~Lineage,scales = "free_x",nrow=4) +
  geom_hline(yintercept = 11) +
  ylab("Variant count 95th percentile")

sublin<-ggplot(data=plot_dist_trop,aes(x=STRAIN,y=p95,color=CHROM,group=CHROM)) + 
  geom_point(size=0.3) + 
  geom_line(linewidth=0.3) + 
  theme_classic() + 
  theme(legend.position = "right",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  facet_wrap(~Sublineage,scales = "free_x",nrow=4) +
  geom_hline(yintercept = 11) +
  ylab("Variant count 95th percentile")

all_percentiles <- rbind(plot_dist,plot_dist_trop %>% dplyr::select(-Sublineage) %>% dplyr::mutate(REF="QX1410")) %>%
  dplyr::select(STRAIN,Lineage) %>%
  dplyr::distinct(STRAIN,.keep_all = T)

excluded <- c("HPT18", "HPT24", "HPT11", "VX34", "JU3326", "NIC893",
              "BRC20503", "BRC20502", "ECA276", "ECA163", "JU2767", "ED3102","QR25")
 
tree_nwk <- ape::read.tree(file="../processed_data/phy_file_LD_0.9.phy.contree")

tree <- ggtree::ggtree(phangorn::midpoint(tree_nwk), layout = "circular")

tip_data <- dplyr::filter(tree$data, isTip) %>%
  dplyr::left_join(all_percentiles, by = c("label" = "STRAIN"))

tree$data <- dplyr::left_join(tree$data, all_percentiles, by = c("label" = "STRAIN")) %>% 
  dplyr::mutate(Lineage=ifelse(label %in% excluded,"Excluded",ifelse(is.na(Lineage),"Reference",Lineage))) %>% 
  dplyr::left_join(lineages %>% dplyr::select(Lineage,lineage_color) %>% dplyr::distinct(Lineage,.keep_all = T),by="Lineage") %>%
  dplyr::mutate(Lineage=ifelse(label=="QX1410"|label=="JU2536","Reference",Lineage)) %>%
  dplyr::mutate(
    lineage_color = dplyr::case_when(
      Lineage == "TH (2)" ~ "#f5f57a",   # lighter version of TH (#aeb400)
      Lineage == "AD (2)" ~ "#ff7a75",   # lighter version of AD (#ff8200)
      Lineage == "Excluded" ~ "#bbbbbb",
      Lineage == "Reference" ~ "#000000",
      TRUE ~ lineage_color)) %>%
  dplyr::mutate(Lineage = factor(Lineage,levels = c("Reference",
                                                    sort(unique(Lineage[!Lineage %in% c("Reference", "Excluded")])),
                                                    "Excluded")))
        
annotree <- tree +
  ggtree::geom_tippoint(ggplot2::aes(color = Lineage), size = 1.5) +
  ggtree::geom_tiplab(
    data = dplyr::filter(tree$data, Lineage == "Reference"),
    ggplot2::aes(label = label),
    size = 4,
    align = TRUE,
    offset = -0.001
  ) +
  theme(legend.position = "right",
        legend.title = element_blank(),
        plot.margin = grid::unit(c(0, 0, 0, 0), "cm")) +
  scale_color_manual(
    values = stats::setNames(tree$data$lineage_color, tree$data$Lineage)
  ) 

ggsave(annotree,filename = "../figures/FigureS13_phylo_annotree_20250930.png", width = 7.5, height = 7.5, device = 'png', dpi = 600, bg = "white")
  
strainL_NR <- unique(SR_stats_WI_NR$STRAIN)

all_stats_NR <- SR_stats_WI_NR %>%
  dplyr::mutate(div_class=ifelse(COUNT >= 11,"I",
                                 ifelse(pc1X < 0.9,"C","R"))) %>%
  dplyr::mutate(div=ifelse(div_class=="C" | div_class=="I","div","nondiv")) %>%               
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::mutate(div_gf=ifelse(div=="nondiv" & lead(div)=="div" & lag(div)=="div","div",div)) %>%
  dplyr::mutate(div_class=ifelse(div==div_gf,div_class,"G")) %>%
  dplyr::ungroup()


# all_stats_NR_freqEst <- frequencyEstimates(all_stats_NR)
regList <- list()
binList <- list()
for (i in 1:length(strainL_NR)) {
  print((i/length(strainL_NR))*100)
  temp <- all_stats_NR %>% dplyr::filter(STRAIN==strainL_NR[i])
  div_call <- clusterBins(temp,"SR")
  div_call[[1]]$STRAIN <- strainL_NR[[i]]
  div_call[[2]]$STRAIN <- strainL_NR[[i]]
  binList[[i]] <- div_call[[1]]
  regList[[i]] <- div_call[[2]]
}

tmp_strains <- temperate_coverage %>% dplyr::select(STRAIN,REF) %>% dplyr::distinct(STRAIN,.keep_all = T)

all_calls_SR_NR <- ldply(regList, data.frame)  %>% dplyr::left_join(lineages,by=c("STRAIN"="isotype"))

gap_clust_NR <- all_calls_SR_NR %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::mutate(forGapSize=lead(minStart)-maxEnd) %>%
  dplyr::mutate(flag3g=ifelse(forGapSize<=5000,"clust","noclust")) %>%
  dplyr::mutate(dec3g=ifelse(flag3g=="clust" ,"join",
                             ifelse(flag3g=="noclust" & lag(flag3g)=="clust","join","nojoin"))) %>%
  dplyr::mutate(dec3g=ifelse(is.na(dec3g),"nojoin",dec3g)) %>%
  dplyr::ungroup()

joinClust_NR <- gap_clust_NR %>% 
  dplyr::filter(dec3g=="join") %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::mutate(segbreak=ifelse(flag3g=="noclust",paste0(dec3g,row_number()),NA)) %>%
  tidyr::fill(segbreak,.direction = 'up') %>%
  dplyr::mutate(gid=data.table::rleid(segbreak)) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(conID=paste0(CHROM,"-",STRAIN,"-",gid)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(conID) %>%
  #dplyr::mutate(gid2=cur_group_id()) #%>%
  dplyr::mutate(newStart=min(minStart),newEnd=max(maxEnd)) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(gapFoot=paste0(rep("R",forGapSize/1e3),collapse = "")) %>%
  dplyr::mutate(new_foot=ifelse(flag3g=="clust",paste0(bin_foot,gapFoot),bin_foot)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(conID) %>%
  dplyr::mutate(clust_foot=paste0(new_foot,collapse = "")) %>%
  dplyr::mutate(newDivSize=newEnd-newStart) %>%
  dplyr::mutate(newMeanVC=mean(meanVC)) %>%
  dplyr::mutate(newMeanCF=mean(meanCF)) %>%
  dplyr::mutate(nclust=n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(conID,.keep_all = T) %>%
  dplyr::select(-minStart,-maxEnd,-divSize,-meanVC,-meanCF,-bin_foot) %>%
  dplyr::rename(minStart=newStart,maxEnd=newEnd,divSize=newDivSize,meanVC=newMeanVC,meanCF=newMeanCF,bin_foot=clust_foot) %>%
  dplyr::select(CHROM,minStart,maxEnd,divSize,meanVC,meanCF,bin_foot,STRAIN,nclust,REF)

nojoin_NR <- gap_clust_NR %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::filter(!(dec3g=="join")) %>%
  dplyr::ungroup() %>%
  dplyr::select(CHROM,minStart,maxEnd,divSize,meanVC,meanCF,bin_foot,STRAIN,REF) %>%
  dplyr::mutate(nclust=1)

all_calls_SR_clustered_NR <- rbind(joinClust_NR,nojoin_NR) %>%
  dplyr::filter(divSize/1e3 >= 5) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ncalls=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(REF) %>%
  dplyr::arrange(REF,desc(ncalls),STRAIN,CHROM,minStart) %>%
  dplyr::mutate(sorter=paste0(ncalls,STRAIN)) %>%
  dplyr::mutate(rleID=data.table::rleid(sorter)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ystrain=cur_group_id()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!grepl("ptg",CHROM))

nuc <- readr::read_tsv("../processed_data/phyloGroup_nucmer_aln2.tsv", col_names = c("S1","E1","S2","E2","L1","L2","IDY","LENR","LENQ","REF","HIFI","STRAIN")) %>%
  dplyr::filter(!grepl("ptg",HIFI)) %>% 
  dplyr::filter(REF==HIFI) %>%
  dplyr::filter(L2 > 3e3) %>%
  dplyr::group_by(STRAIN,HIFI) %>%
  dplyr::arrange(STRAIN,HIFI,S2) %>%
  dplyr::mutate(leadS1=lead(S1),leadE1=lead(E1),leadS2=lead(S2),leadE2=lead(E2),lagS1=lag(S1),lagE1=lag(E1),lagS2=lag(S2),lagE2=lag(E2)) %>%
  dplyr::ungroup()

############# DT test
calls <- as.data.table(all_calls_SR_clustered_NR)[, rowid := .I]
calls <- calls[, .(rowid, region_source = STRAIN, genome = REF, contig = CHROM,
                   start = minStart, end = maxEnd)]

nuc_dt <- as.data.table(nuc)[, idx := .I]
nuc_dt <- nuc_dt[, .(idx, genome = STRAIN, contig = HIFI,
                     start = pmin(S2, E2), end = pmax(S2, E2),
                     orig_start=S2,orig_end=E2, L1=L1, L2=L2,
                     leadS1=leadS1,leadE1=leadE1,leadS2=leadS2,leadE2=leadE2,
                     lagS1=lagS1,lagE1=lagE1,lagS2=lagS2,lagE2=lagE2,
                     refchrom=REF,refstart=S1,refend=E1)]

# Remove old values in case you're re-running in the same session
#nuc_dt[, `:=`(region_id = NULL, region_source = NULL)]
calls[, `:=`(start = as.integer(start), end = as.integer(end))]
nuc_dt[, `:=`(start = as.integer(start), end = as.integer(end))]
# Set proper keys
setkey(calls, genome, contig, start, end)
setkey(nuc_dt, genome, contig, start, end)

# Rerun safe and clean overlap
matched <- foverlaps(
  x = calls,
  y = nuc_dt,
  type = "any",
  nomatch = 0
)

matched_df <- as.data.frame(matched) %>%
  dplyr::mutate(hdr_chrom=contig) %>%
  dplyr::mutate(INV=ifelse(orig_start==start,F,T)) %>%
  dplyr::select(refstart,refend,orig_start,orig_end,L1,L2,refchrom,contig,genome,INV,start,end,rowid,region_source,hdr_chrom,i.start,i.end,leadS1,leadE1,leadS2,leadE2,lagS1,lagE1,lagS2,lagE2) %>%
  dplyr::rename(S1=refstart,E1=refend,S2=orig_start,E2=orig_end,REF=refchrom,HIFI=contig,HIFI_strain=genome,St2=start,Et2=end,group_id=rowid,hdr_strain=region_source,hdr_start=i.start,hdr_end=i.end) %>% 
  dplyr::arrange(group_id,S2) %>%
  dplyr::mutate(HDRid = paste0(hdr_strain,hdr_chrom,hdr_start,hdr_end))


tigFilt2 <- matched_df %>%
  dplyr::arrange(group_id,S1) %>%
  dplyr::group_by(group_id) %>%
  dplyr::mutate(leadDiff=lead(S1)-E1) %>%
  dplyr::mutate(jump=ifelse(leadDiff > 5E4,1,0)) %>%
  dplyr::mutate(leadDiff=ifelse(is.na(leadDiff),0,leadDiff)) %>%
  dplyr::mutate(run_id = cumsum(c(1, head(jump, -1)))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(group_id,run_id) %>%
  dplyr::mutate(gsize=n()) %>%
  dplyr::mutate(len=abs(E1-S1)) %>%
  dplyr::mutate(sumlen=sum(len)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(group_id) %>%
  dplyr::filter(sumlen==max(sumlen)) %>%
  dplyr::select(-gsize) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group_id,S1) 


# ggplot(tigFilt2 %>% dplyr::filter(HDRid == "TWN1824X932000942000")) +
#   geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=hdr_start,ymax=hdr_end),fill="lightgrey") +
#   geom_segment(aes(x=S1,xend=E1,y=S2,yend=E2)) 

trim_spacer = 1e3
#trims long alignments to the focal region (i.e. hap_start to hap_end, but transformed to the other genome)
tigTrim <- tigFilt2 %>%
  dplyr::group_by(group_id) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(scale_distortion = ((L2 - L1)/L1)) %>%
  dplyr::mutate(rboundDist=hdr_start-min(S2,E2)) %>%
  #dplyr::mutate(E1=ifelse(rboundDist>trim_spacer & INV==F,(E1-(rboundDist-trim_spacer)),E1)) %>%
  dplyr::mutate(S1=ifelse(rboundDist>trim_spacer & INV==F ,(S1+(rboundDist-trim_spacer+(rboundDist*scale_distortion))),S1)) %>%
  dplyr::mutate(E1=ifelse(rboundDist>trim_spacer & INV==T ,(E1-(rboundDist-trim_spacer+(rboundDist*scale_distortion))),E1)) %>%
  dplyr::mutate(S2=ifelse(rboundDist>trim_spacer & INV==F,(S2+(rboundDist-trim_spacer)),S2)) %>%
  dplyr::mutate(E2=ifelse(rboundDist>trim_spacer & INV==T,(E2+(rboundDist-trim_spacer)),E2)) %>%
  dplyr::mutate(lboundDist=max(S2,E2)-hdr_end) %>%
  dplyr::mutate(E1=ifelse(lboundDist>trim_spacer & INV==F ,(E1-(lboundDist-trim_spacer+(lboundDist*scale_distortion))),E1)) %>%
  dplyr::mutate(S1=ifelse(lboundDist>trim_spacer & INV==T ,(S1+(lboundDist-trim_spacer+(lboundDist*scale_distortion))),S1)) %>%
  dplyr::mutate(E2=ifelse(lboundDist>trim_spacer & INV==F,(E2-(lboundDist-trim_spacer)),E2)) %>%
  dplyr::mutate(S2=ifelse(lboundDist>trim_spacer & INV==T,(S2-(lboundDist-trim_spacer)),S2)) %>%
  dplyr::ungroup()
  
# ggplot(tigTrim %>% dplyr::filter(HDRid == "TWN1824X932000942000")) +
#   geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=hdr_start,ymax=hdr_end),fill="lightgrey") +
#   #geom_segment(aes(x=leadS1,xend=leadE1,y=leadS2,yend=leadE2, color="outreg_lead")) +
#   #geom_segment(aes(x=lagS1,xend=lagE1,y=lagS2,yend=lagE2, color="outreg_lag")) +
#   geom_segment(aes(x=S1,xend=E1,y=S2,yend=E2)) +
#   xlab("REF coords") +
#   ylab("WI coords")


# ggplot(tigTrim %>% dplyr::filter(HDRid == "TWN1824X932000942000")) +
#   geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=hdr_start,ymax=hdr_end),fill="lightgrey") +
#   geom_segment(aes(x=leadS1,xend=leadE1,y=leadS2,yend=leadE2, color="outreg_lead")) +
#   geom_segment(aes(x=lagS1,xend=lagE1,y=lagS2,yend=lagE2, color="outreg_lag")) +
#   geom_segment(aes(x=S1,xend=E1,y=S2,yend=E2, color="inreg_trimmed")) +
#   xlab("REF coords") +
#   ylab("WI coords") 

tigMarkExtend <- tigTrim %>%
  dplyr::group_by(group_id) %>%
  dplyr::mutate(E2_extend=ifelse(INV==F & E2 == max(E2) & E2 < hdr_end, T, F),
                S2_extend=ifelse(INV==F & S2 == min(S2) & S2 > hdr_start, T,F),
                iE2_extend=ifelse(INV==T & E2 == min(E2) & E2 > hdr_start,T,F),
                iS2_extend=ifelse(INV==T & S2 == max(S2) & S2 < hdr_end, T,F)) %>%
  dplyr::mutate(any_extend=ifelse(E2_extend == T | S2_extend == T | iE2_extend==T | iS2_extend ==T,T,F)) %>%
  dplyr::ungroup()
 

tigToExtend <- tigMarkExtend %>% 
  dplyr::filter(any_extend==T) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(extend_length_WI_lead=ifelse(E2_extend==T & min(leadS2,leadE2) > hdr_end, min(leadS2,leadE2)-E2,
                                             ifelse(iS2_extend==T & min(leadS2,leadE2) > hdr_end, min(leadS2,leadE2)-S2,NA)),
                extend_length_WI_lag=ifelse(S2_extend==T & max(lagS2,lagE2) < hdr_start, S2-max(lagS2,lagE2),
                                               ifelse(iE2_extend==T & max(lagS2,lagE2) < hdr_start, E2-max(lagS2,lagE2),NA)),
                extend_length_REF_lead=ifelse(E2_extend==T & min(leadS2,leadE2) > hdr_end, 
                                              ifelse(leadS1 >= E1,leadS1-E1,ifelse(leadE1>=E1,leadE1-E1,NA)),
                                              ifelse(iS2_extend==T & min(leadS2,leadE2) > hdr_end, ifelse(S1>=leadE1,S1-leadE1,ifelse(S1>=leadS1,S1-leadS1,NA)),NA)),
                extend_length_REF_lag=ifelse(S2_extend==T & max(lagS2,lagE2) < hdr_start, 
                                            ifelse(lagE1<=S1,S1-lagE1,ifelse(lagS1<=S1,S1-lagS1,NA)),
                                            ifelse(iE2_extend==T & max(lagS2,lagE2) < hdr_start, ifelse(lagS1>=E1,lagS1-E1,ifelse(lagE1>=E1,lagE1-E1,NA)),NA))) %>%
  dplyr::ungroup()


extendDat <- rbind(tigToExtend %>% 
                     dplyr::select(extend_length_REF_lead,extend_length_WI_lead) %>% 
                     dplyr::rename(extend_length_WI=extend_length_WI_lead,extend_length_REF=extend_length_REF_lead),
                   tigToExtend %>% 
                     dplyr::select(extend_length_REF_lag,extend_length_WI_lag) %>% 
                     dplyr::rename(extend_length_WI=extend_length_WI_lag,extend_length_REF=extend_length_REF_lag)) %>%
  dplyr::filter(!is.na(extend_length_WI) & !is.na(extend_length_REF))

sc <- ggplot(data=extendDat) + 
  geom_point(aes(x=extend_length_REF/1e3,y=extend_length_WI/1e3),size=1) + 
  geom_rect(xmin=0,xmax=100,ymin=-Inf,ymax=Inf,fill=NA,color="grey",linetype="dashed")+
  geom_rect(xmin=Inf,xmax=-Inf,ymin=0,ymax=100,fill=NA,color="grey",linetype="dashed")+
  theme_classic() + 
  xlab("REF extension distances (kb)") + 
  ylab("WI extension distances (kb)") +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_x_continuous(expand = c(0.01,0))

h1 <- ggplot(data=extendDat) + 
  geom_histogram(aes(x=extend_length_REF/1e3),binwidth = 1) + 
  theme_classic() + 
  xlab("") + 
  ylab("count") + 
  coord_cartesian(xlim=c(0,100)) +
  scale_y_continuous(expand = c(0.01,0),
                     labels = function(y) y / 1000,
                     name = "count (thousand)") +
  scale_x_continuous(expand = c(0.01,0))

h2 <- ggplot(data=extendDat) + 
  geom_histogram(aes(y=extend_length_WI/1e3),binwidth = 1) + 
  theme_classic() + ylab("") + 
  coord_cartesian(ylim=c(0,100)) +
  scale_x_continuous(labels = function(x) x / 1000,
                     expand = c(0.01, 0),
                     name = "count (thousand)") +
  scale_y_continuous(expand = c(0.01,0))


middle_row <- cowplot::plot_grid(
  sc,       
  h2, 
  ncol = 2,
  rel_widths = c(4, 1), 
  align = "hv"
)

empty_plot <- ggplot() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_void() 

top_row <- plot_grid(
  empty_plot,
  h1,
  empty_plot,
  ncol = 3,
  rel_widths = c(0.09,4, 1),
  align = "hv"
)

final_plot <- cowplot::plot_grid(
  top_row,      
  middle_row,    
  ncol = 1,
  rel_heights = c(1, 4),  
  align = "v"
)

tigExtensions <- rbind(tigToExtend,tigMarkExtend %>% dplyr::filter(any_extend==F) %>% dplyr::mutate(extend_length_WI_lead=NA,extend_length_REF_lead=NA,extend_length_WI_lag=NA,extend_length_REF_lag=NA))

ggplot(tigExtensions %>% dplyr::filter(HDRid == "TWN1824X932000942000")) +
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=hdr_start,ymax=hdr_end),fill="lightgrey") +
  geom_segment(aes(x=leadS1,xend=leadE1,y=leadS2,yend=leadE2, color="outreg_lead")) +
  geom_segment(aes(x=lagS1,xend=lagE1,y=lagS2,yend=lagE2, color="outreg_lag")) +
  geom_segment(aes(x=S1,xend=E1,y=S2,yend=E2, color="inreg_trimmed")) +
  xlab("REF coords") +
  ylab("WI coords")

tigExtended_50kb <- tigExtensions %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(E1=ifelse(E2_extend==T & min(leadS2,leadE2) > hdr_end & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < 5e4 & extend_length_REF_lead < 5e4, ifelse(leadS1 >= E1,leadS1,ifelse(leadE1>=E1,leadE1,E1)),E1)) %>%
  dplyr::mutate(E1=ifelse(iE2_extend==T & max(lagS2,lagE2) < hdr_start & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < 5e4 & extend_length_REF_lag < 5e4, ifelse(lagS1>=E1,lagS1,ifelse(lagE1>=E1,lagE1,E1)),E1)) %>%
  dplyr::mutate(S1=ifelse(S2_extend==T & max(lagS2,lagE2) < hdr_start & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < 5e4 & extend_length_REF_lag < 5e4, ifelse(lagE1<=S1,lagE1, ifelse(lagS1<=S1,lagS1,S1)),S1)) %>%
  dplyr::mutate(S1=ifelse(iS2_extend==T & min(leadS2,leadE2) > hdr_end & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < 5e4 & extend_length_REF_lead < 5e4, ifelse(S1>=leadE1,leadE1,ifelse(S1>=leadS1,leadS1,S1)),S1)) %>%
  dplyr::mutate(E2=ifelse(E2_extend==T & min(leadS2,leadE2) > hdr_end & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < 5e4 & extend_length_REF_lead < 5e4, min(leadS2,leadE2),E2)) %>%
  dplyr::mutate(E2=ifelse(iE2_extend==T & max(lagS2,lagE2) < hdr_start & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < 5e4 & extend_length_REF_lag < 5e4, max(lagS2,lagE2),E2)) %>%
  dplyr::mutate(S2=ifelse(S2_extend==T & max(lagS2,lagE2) < hdr_start & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < 5e4 & extend_length_REF_lag < 5e4, max(lagS2,lagE2),S2)) %>%
  dplyr::mutate(S2=ifelse(iS2_extend==T & min(leadS2,leadE2) > hdr_end & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < 5e4 & extend_length_REF_lead < 5e4, min(leadS2,leadE2),S2)) %>% 
  dplyr::ungroup()

tigExtended_25kb <- tigExtensions %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(E1=ifelse(E2_extend==T & min(leadS2,leadE2) > hdr_end & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < 2.5e4 & extend_length_REF_lead < 2.5e4, ifelse(leadS1 >= E1,leadS1,ifelse(leadE1>=E1,leadE1,E1)),E1)) %>%
  dplyr::mutate(E1=ifelse(iE2_extend==T & max(lagS2,lagE2) < hdr_start & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < 2.5e4 & extend_length_REF_lag < 2.5e4, ifelse(lagS1>=E1,lagS1,ifelse(lagE1>=E1,lagE1,E1)),E1)) %>%
  dplyr::mutate(S1=ifelse(S2_extend==T & max(lagS2,lagE2) < hdr_start & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < 2.5e4 & extend_length_REF_lag < 2.5e4, ifelse(lagE1<=S1,lagE1, ifelse(lagS1<=S1,lagS1,S1)),S1)) %>%
  dplyr::mutate(S1=ifelse(iS2_extend==T & min(leadS2,leadE2) > hdr_end & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < 2.5e4 & extend_length_REF_lead < 2.5e4, ifelse(S1>=leadE1,leadE1,ifelse(S1>=leadS1,leadS1,S1)),S1)) %>%
  dplyr::mutate(E2=ifelse(E2_extend==T & min(leadS2,leadE2) > hdr_end & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < 2.5e4 & extend_length_REF_lead < 2.5e4, min(leadS2,leadE2),E2)) %>%
  dplyr::mutate(E2=ifelse(iE2_extend==T & max(lagS2,lagE2) < hdr_start & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < 2.5e4 & extend_length_REF_lag < 2.5e4, max(lagS2,lagE2),E2)) %>%
  dplyr::mutate(S2=ifelse(S2_extend==T & max(lagS2,lagE2) < hdr_start & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < 2.5e4 & extend_length_REF_lag < 2.5e4, max(lagS2,lagE2),S2)) %>%
  dplyr::mutate(S2=ifelse(iS2_extend==T & min(leadS2,leadE2) > hdr_end & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < 2.5e4 & extend_length_REF_lead < 2.5e4, min(leadS2,leadE2),S2)) %>% 
  dplyr::ungroup()

# id="ECA287II1533900015346000"
# id="QG2996V1775500017792000"
# id="NIC814V1813100018140000"
# id="NIC1631X75060007539000"
# id="QG1126V1696300017079000"
# id="QG4093II847000855000"
# ggplot(rbind(tigExtensions %>% dplyr::filter(HDRid == id) %>% dplyr::mutate(type="original"),tigExtended_25kb %>% dplyr::filter(HDRid == id) %>% dplyr::mutate(type="extended_max25kb"),tigExtended_50kb %>% dplyr::filter(HDRid == id) %>% dplyr::mutate(type="extended_max50kb"))) +
#   geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=hdr_start,ymax=hdr_end),fill="lightgrey") +
#   geom_segment(aes(x=leadS1,xend=leadE1,y=leadS2,yend=leadE2, color="outreg_lead")) +
#   geom_segment(aes(x=lagS1,xend=lagE1,y=lagS2,yend=lagE2, color="outreg_lag")) +
#   geom_segment(aes(x=S1,xend=E1,y=S2,yend=E2, color="inreg_trimmed")) +
#   xlab("REF coords") +
#   ylab("WI coords") +
#   facet_wrap(~as.factor(type),scales = 'free')

counts25kb <- tigExtended_25kb %>%
  filter(
    any_extend == TRUE,
    (!is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag)) | (!is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead)),
    (extend_length_WI_lag < 2.5e4 & extend_length_REF_lag < 2.5e4) | (extend_length_WI_lead < 2.5e4 & extend_length_REF_lead < 2.5e4)
  ) %>%
  group_by(hdr_strain, HDRid) %>%
  summarise(has_extension = any(any_extend), .groups = "drop") %>%  # optional, since filter already ensures this
  group_by(hdr_strain) %>%
  summarise(count_true = n(), .groups = "drop") %>%
  mutate(window_size = "50kb")

counts50kb <- tigExtended_50kb %>%
  filter(
    any_extend == TRUE,
    (!is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag)) | (!is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead)),
    (extend_length_WI_lag < 5e4 & extend_length_REF_lag < 5e4) | (extend_length_WI_lead < 5e4 & extend_length_REF_lead < 5e4)
  ) %>%
  group_by(hdr_strain, HDRid) %>%
  summarise(has_extension = any(any_extend), .groups = "drop") %>%  # optional, since filter already ensures this
  group_by(hdr_strain) %>%
  summarise(count_true = n(), .groups = "drop") %>%
  mutate(window_size = "50kb")

hdr_counts <- tigTrim %>%
  group_by(hdr_strain) %>%
  summarise(num_unique_HDRid = n_distinct(HDRid), .groups = "drop")

combined_counts <- counts25kb %>% dplyr::left_join(counts50kb,by=c("hdr_strain"))  %>% dplyr::left_join(hdr_counts,by=c("hdr_strain"))%>%
  dplyr::arrange(num_unique_HDRid) %>%
  dplyr::mutate(hdr_strain = factor(hdr_strain, levels = hdr_strain)) %>%
  dplyr::mutate(diff=count_true.y-count_true.x)

# ggplot(combined_counts) +
#   geom_bar(aes(x = hdr_strain, y = num_unique_HDRid, fill="all_HDRs"),color="black",stat = "identity",position = position_identity()) +
#   geom_bar(aes(x = hdr_strain, y = count_true.y, fill="HDRs_ext_50kb"),color="black",stat = "identity",position = position_identity()) +
#   geom_bar(aes(x = hdr_strain, y = count_true.x, fill="HDRs_ext_25kb"),color="black",stat = "identity",position = position_identity()) +  # or: position = "stack"
#   labs(x = "Strain", y = "Count of HDRs", fill = "") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#   scale_fill_manual(values = c("all_HDRs"="grey","HDRs_ext_50kb"="blue","HDRs_ext_25kb"="forestgreen"))

hdr_transformed_orig <- tigTrim %>%
  dplyr::group_by(group_id) %>%
  dplyr::summarise(
    S1 = min(S1, E1),
    E1 = max(S1, E1),
    across(
      .cols = -c(S1, E1, S2, E2, St2, Et2),
      .fns = dplyr::first
    ),
    .groups = "drop"
  )

 # hdr_transformed_25ext <- tigExtended_25kb %>%
 #    group_by(group_id) %>%
 #    summarise(
 #      S1 = min(S1, E1),
 #      E1 = max(S1, E1),
 #      across(
 #        .cols = -c(S1, E1, S2, E2, St2, Et2),
 #        .fns = first
 #      ),
 #      .groups = "drop"
 #    )
 
 hdr_transformed_50ext <- tigExtended_50kb %>%
   dplyr::group_by(group_id) %>%
   dplyr::summarise(
     S1 = min(S1, E1),
     E1 = max(S1, E1),
     across(
       .cols = -c(S1, E1, S2, E2, St2, Et2),
       .fns = dplyr::first
     ),
     .groups = "drop"
   )

 
 gap_clust_NR_TR <- hdr_transformed_50ext %>%
   dplyr::select(REF,S1,E1,hdr_strain,HIFI_strain) %>%
   dplyr::rename(CHROM=REF,minStart=S1,maxEnd=E1,STRAIN=hdr_strain,REF=HIFI_strain) %>%
   dplyr::mutate(divSize=maxEnd-minStart) %>%
   dplyr::arrange(STRAIN,CHROM,minStart) %>%
   dplyr::group_by(STRAIN,CHROM) %>%
   dplyr::mutate(forGapSize=lead(minStart)-maxEnd) %>%
   dplyr::mutate(flag3g=ifelse(forGapSize<=5000,"clust","noclust")) %>%
   dplyr::mutate(dec3g=ifelse(flag3g=="clust" ,"join",
                              ifelse(flag3g=="noclust" & lag(flag3g)=="clust","join","nojoin"))) %>%
   dplyr::mutate(dec3g=ifelse(is.na(dec3g),"nojoin",dec3g)) %>%
   dplyr::ungroup()
 
 joinClust_NR_TR <- gap_clust_NR_TR %>% 
   dplyr::filter(dec3g=="join") %>%
   dplyr::group_by(STRAIN,CHROM) %>%
   dplyr::mutate(segbreak=ifelse(flag3g=="noclust",paste0(dec3g,row_number()),NA)) %>%
   tidyr::fill(segbreak,.direction = 'up') %>%
   dplyr::mutate(gid=data.table::rleid(segbreak)) %>%
   dplyr::ungroup() %>%
   dplyr::rowwise() %>%
   dplyr::mutate(conID=paste0(CHROM,"-",STRAIN,"-",gid)) %>%
   dplyr::ungroup() %>%
   dplyr::group_by(conID) %>%
   dplyr::mutate(newStart=min(minStart),newEnd=max(maxEnd)) %>%
   dplyr::ungroup() %>%
   dplyr::group_by(conID) %>%
   dplyr::mutate(newDivSize=newEnd-newStart) %>%
   dplyr::mutate(nclust=n()) %>%
   dplyr::ungroup() %>%
   dplyr::distinct(conID,.keep_all = T) %>%
   dplyr::select(-minStart,-maxEnd,-divSize) %>%
   dplyr::rename(minStart=newStart,maxEnd=newEnd,divSize=newDivSize) %>%
   dplyr::select(CHROM,minStart,maxEnd,divSize,STRAIN,nclust,REF)
 
 nojoin_NR_TR <- gap_clust_NR_TR %>%
   dplyr::group_by(STRAIN,CHROM) %>%
   dplyr::filter(!(dec3g=="join")) %>%
   dplyr::ungroup() %>%
   dplyr::select(CHROM,minStart,maxEnd,divSize,STRAIN,REF) %>%
   dplyr::mutate(nclust=1)
 
 all_calls_SR_clustered_NR_TR <- rbind(joinClust_NR_TR,nojoin_NR_TR) %>%
   dplyr::filter(divSize/1e3 >= 5) %>%
   dplyr::group_by(STRAIN) %>%
   dplyr::mutate(ncalls=n()) %>%
   dplyr::ungroup() %>%
   dplyr::group_by(REF) %>%
   dplyr::arrange(REF,desc(ncalls),STRAIN,CHROM,minStart) %>%
   dplyr::mutate(sorter=paste0(ncalls,STRAIN)) %>%
   dplyr::mutate(rleID=data.table::rleid(sorter)) %>%
   dplyr::ungroup() %>%
   dplyr::group_by(STRAIN) %>%
   dplyr::mutate(ystrain=cur_group_id()) %>%
   dplyr::ungroup() %>%
   dplyr::filter(!grepl("ptg",CHROM))

hdr_rest_unt <- all_calls_SR_clustered_NR %>%
  dplyr::select(CHROM,minStart,maxEnd,STRAIN,REF) %>%
  dplyr::rename(source=REF) %>%
  dplyr::mutate(mode="UNTRANSFORMED")

hdr_rest <- all_calls_SR_clustered_NR_TR %>%
  dplyr::select(CHROM,minStart,maxEnd,STRAIN,REF) %>%
  dplyr::rename(source=REF) %>%
  dplyr::mutate(mode="TRANSFORMED")


hdr_rest_comp <- rbind(hdr_rest_unt,hdr_rest) %>%
  dplyr::mutate(divSize=maxEnd-minStart) %>%
  dplyr::filter(divSize >= 5e3) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ncalls=n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(mode,desc(source),desc(ncalls),STRAIN,CHROM,minStart) %>%
  dplyr::group_by(mode) %>%
  dplyr::mutate(sorter=paste0(ncalls,STRAIN)) %>%
  dplyr::mutate(rleID=data.table::rleid(sorter)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ystrain=cur_group_id()) %>%
  dplyr::ungroup()

restcomp <- ggplot(hdr_rest_comp %>% dplyr::filter(divSize >= 5e3)) + 
  geom_rect(aes(xmin=minStart/1e6,xmax=maxEnd/1e6,ymin=rleID-0.45,ymax=rleID+0.45,fill=source)) + 
  facet_grid(mode~CHROM, scales = 'free_x') + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        panel.border = element_rect(fill = NA),
        axis.ticks.y = element_blank(),
        axis.text.y=element_blank())+
  #strip.text.x = element_blank())  +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Physical position (Mb)") +
  ylab("280 Isotype strains")

hdr_qx <- all_calls_SR_clustered %>%
  dplyr::select(CHROM,minStart,maxEnd,STRAIN) %>%
  dplyr::mutate(source="QX1410")

hdr_tot <- rbind(hdr_qx,hdr_rest %>% dplyr::select(-mode) %>% dplyr::mutate(minStart=round(minStart),maxEnd=round(maxEnd))) %>%
  dplyr::mutate(divSize=maxEnd-minStart) %>%
  dplyr::filter(divSize >= 5e3) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ncalls=n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(ncalls),STRAIN,CHROM,minStart) %>%
  dplyr::mutate(sorter=paste0(ncalls,STRAIN)) %>%
  dplyr::mutate(rleID=data.table::rleid(sorter)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ystrain=cur_group_id()) %>%
  dplyr::ungroup()

options(scipen=10)
write.table(hdr_tot %>%
              dplyr::select(CHROM,minStart,maxEnd,STRAIN,source,divSize) %>%
              dplyr::arrange(STRAIN,CHROM,minStart,maxEnd),
            file="../tables/TableS7_HDR_allStrain_5kbclust_20250930.tsv",row.names = F,quote = F,sep = '\t')
write.table(hdr_rest_unt %>%
              dplyr::mutate(divSize=maxEnd-minStart) %>%
              dplyr::select(CHROM,minStart,maxEnd,STRAIN,source,divSize) %>%
              dplyr::arrange(STRAIN,CHROM,minStart,maxEnd),
            file="../processed_data/HDR_CB_otherRG_UNT_5kbclust_20250930.tsv",row.names = F,quote = F,sep = '\t')
options(scipen=0)


