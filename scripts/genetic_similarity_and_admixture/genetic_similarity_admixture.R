library(dplyr)
library(tidyr)
library(readr)
library(ape)
library(viridis)
library(pheatmap)
library(gplots)
library(gridExtra)
library(grid)
library(ggplot2)
library(cowplot)
library(ggplotify)
library(phytools)
library(geosphere)
library(igraph)
library(ggraph)
library(ggforce)
library(ComplexHeatmap)
library(circlize)
library(admixtools)

#set geographic region colors
geo_colors <- c("Hawaii"="#66C2A5", 
                "Australia"="#FC8D62", 
                "Central America"="#8DA0CB",
                "South America"="#E78AC3", 
                "Africa"="#A6D854", 
                "Caribbean"="#FFD92F",
                "Taiwan" = "#E5C494", 
                "North America" = "#A65628", 
                "Europe" = "#E41A1C",
                "Asia" = "#399EB8", 
                "Pacific" = "#611EA1",
                "New Zealand" = "green",
                "Atlantic" = "purple", 
                "Oceania" ="#DB7779", 
                "Micronesia" = "#E7298A",
                "Indonesia" = "#7570B3", 
                "Malay Archipelago" = "#4110B3", 
                "C. America"="#8DA0CB",
                "S. America"="#E78AC3", 
                "N. America" = "#A65628", 
                "Cosmopolitan" = "gray30",
                "unknown" = 'grey')

#make geo colors data frame
df_colors <- data.frame(unname(geo_colors),names(geo_colors)) %>% dplyr::rename(color=`unname.geo_colors.`, geo=`names.geo_colors.`)

lineages <- readr::read_tsv("../../processed_data/genetic_similarity_and_admixutre/isotype_byRG_GeoLocAdmCol_20250909.tsv") %>%
  dplyr::mutate(sublineage_color=ifelse(Sublineage=="TC","#ff0000",sublineage_color)) %>%
  dplyr::mutate(Sublineage=ifelse(Sublineage=="TC","TT",Sublineage))  %>%
  dplyr::filter(Lineage!="FM")

#read isotype assignments
isos <- readr::read_tsv(file="../../processed_data/genetic_similarity_and_admixutre/isotype_groups.tsv") %>%
  dplyr::group_by(isotype) %>%
  dplyr::summarise(count=n())

geo_info <- readr::read_csv(file="../../processed_data/genetic_similarity_and_admixutre/Cb_indep_isotype_info_geo.csv") %>%
  dplyr::left_join(df_colors,by=c("geo")) %>%
  dplyr::left_join(lineages %>% dplyr::select(isotype, Lineage, Sublineage,lineage_color,sublineage_color), by="isotype") %>%
  dplyr::filter(!is.na(Lineage)) %>%
  dplyr::mutate(abslat=abs(lat)) 

cvmat <- readr::read_tsv(file="../../processed_data/genetic_similarity_and_admixutre/cv_matrix.tsv")

cvmat_long <- cvmat %>%
  tidyr::pivot_longer(
    cols = starts_with("S"),      # all seed columns
    names_to = "seed",
    values_to = "cv_error",
    values_drop_na = T) %>%
  dplyr::mutate(seed = as.integer(sub("^S", "", seed)))  # strip "S" and make numeric

### Any K19+ seems viable
cv <- ggplot() +
  geom_jitter(data=cvmat_long,aes(x=K,y=cv_error,group=K),alpha=0.7,size=0.5, width = 0.2, height = 0) +
  geom_boxplot(data=cvmat_long,aes(x=K,y=cv_error,group=K),fill=NA,outliers = F) +
  #coord_cartesian(ylim=c(0,0.1))+
  labs(color="Seed") +
  scale_x_continuous(breaks = seq(min(cvmat_long$K), max(cvmat_long$K), by = 1), expand = c(0,0)) +
  ylim(0,NA)+
  theme_bw() +
  ylab("CV")
  
#lets compare K2-30 to our relatedness groups
kmin=2
kmax=30
fraclist <- list()
qlist <- list()
for (i in kmin:kmax) {
    #if reach maximum of LETTERS[], then increase to two-letter subpop names
    if (i>26) {
      colnames = c()
      for (i in 27:i){
        symbol=paste0(LETTERS[floor((i+1)/26)],LETTERS[i-(26*floor((i+1)/26))])
        colnames=c(colnames,symbol)
      }
      colnames=c(LETTERS[1:26],colnames)
    } else {
      colnames=LETTERS[1:i]
      symbol=LETTERS[i]
    }
    
    tmpQ <- readr::read_tsv(paste0("../../processed_data/genetic_similarity_and_admixutre/concat_Qfiles_K",i,".tsv"),col_names=c("isotype",colnames,"run_ID"))
  
    qlist[[i-(kmin-1)]] <- tmpQ %>%
      tidyr::pivot_longer(
        cols = A:!!sym(symbol),
        names_to = "subpop",
        values_to = "fraction") %>% 
      tidyr::separate(run_ID, into=c("LD_tag","LD_val","K_val","seed"),sep = "_") %>%
      dplyr::select(-LD_tag) %>%
      dplyr::left_join(geo_info,by="isotype")
    
  fraclist[[i-(kmin-1)]] <- qlist[[i-(kmin-1)]] %>%
    dplyr::filter(fraction >= 0.999) #get non-admixed individuals
}

ctlist<- list()
rgm_perseed <- list()
for (i in kmin:kmax) {
  RGassignment<-fraclist[[i-(kmin-1)]] %>%
    dplyr::group_by(K_val,seed) %>%
    dplyr::mutate(nlin_assigned=length(unique(Lineage))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(K_val,seed,subpop) %>%
    dplyr::mutate(nlin=length(unique(Lineage))) %>%
    dplyr::distinct(Lineage,.keep_all = T) %>%
    dplyr::select(-isotype) %>%
    dplyr::ungroup() 
  
  mergedRG<- RGassignment %>%  dplyr::filter(nlin>1) 
  
  ctlist[[i-(kmin-1)]] <- c(i,length(unique(mergedRG$seed)),max(mergedRG$nlin_assigned),max(mergedRG$nlin))
  
   per_seed <- RGassignment %>% 
     dplyr::group_by(K_val,seed) %>%
     dplyr::filter(nlin==max(nlin)) %>%
     dplyr::distinct(seed,.keep_all = T)
  
  rgm_perseed[[i-(kmin-1)]] <- per_seed$nlin
}
df <- as.data.frame(do.call(rbind, ctlist)) %>% dplyr::filter(V1>=10)
df_perseed <- do.call(rbind, lapply(seq_along(rgm_perseed), function(i) {data.frame(K   = i + (kmin - 1), num = rgm_perseed[[i]])})) %>% dplyr::filter(K>=10)
colnames(df) <- c("K", "nrun_with_merged_RG", "max_nonadmix_RGs", "max_nonadmix_RG_merged")

rg1 <- ggplot()+geom_point(data=df,aes(x=K,y=nrun_with_merged_RG),size=0.5) + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5),
        axis.text = element_text(size=6),
        axis.title = element_text(size=9),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size=7)) +
  ylab("Number of independent runs\nwith multiple RGs per subpopulation") + 
  scale_x_continuous(breaks = seq(min(df$K), max(df$K), by = 1)) + 
  scale_y_continuous(breaks = seq(min(df$nrun_with_merged_RG), max(df$nrun_with_merged_RG), by = 1)) 

rg2 <- ggplot()+geom_point(data=df,aes(x=K,y=max_nonadmix_RGs),size=0.5) + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5),
        axis.text = element_text(size=6),
        axis.title = element_text(size=9),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size=7)) +
  ylab("Maximum number of  RGs\nwith non-admixed individuals") + 
  scale_x_continuous(breaks = seq(min(df$K), max(df$K), by = 1)) + 
  scale_y_continuous(breaks = seq(min(df$max_nonadmix_RGs), max(df$max_nonadmix_RGs), by = 1))

df_counts <- df_perseed %>%
  dplyr::count(K, num)

rg3<-ggplot(df_counts, aes(x = K, y = num, fill = n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), color = "black",size=2.5) +
  scale_fill_gradient(low = "lightblue", high = "blue") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5),
        axis.text = element_text(size=6),
        axis.title = element_text(size=9),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size=7),
        legend.justification.right = "top",
        legend.margin = margin(0,0,0,0),
        plot.margin = margin(0,0,0,0))+
  labs(fill = "Number of\nindependent\nruns") +   # <-- this sets legend title
  ylab("Maximum number of RGs\nper subpopulation") +
  scale_x_continuous(breaks = seq(min(df_perseed$K), max(df_perseed$K), by = 1)) +
  scale_y_continuous(breaks = seq(min(df_perseed$num), max(df_perseed$num), by = 1)) 
  
comprg <- cowplot::plot_grid(rg2,rg1,rg3,nrow=1,ncol=3,rel_widths = c(0.9,0.9,1.2), align="h",axis = "tb", labels=c("b","c","d"))
rg_cv <- cowplot::plot_grid(cv + theme(panel.grid.major = element_line(color="grey80"),panel.grid.minor = element_blank()),comprg,nrow=2,rel_heights = c(1.2,1), align = "v",axis="r",labels=c("a",NA))
ggsave(rg_cv,file="../../figures/SF4_ADX_rg_cv.png",width = 7,height = 6,units = "in",device = 'png',bg="white",dpi=600)

best_kval = 22 #  K=22 selected from plot above (highest assignment stability)
best_k <- cvmat_long %>% dplyr::filter(K==best_kval) %>% dplyr::filter(cv_error==min(cv_error)) #seed 1553 has min cv_error,
best_NADM_assignments <- fraclist[[best_kval-1]] %>% dplyr::filter(seed==(best_k %>% dplyr::pull(seed))) 
best_qfile <- qlist[[best_kval-1]] %>% dplyr::filter(seed==(best_k %>% dplyr::pull(seed)))

#write cluster fractions of selected K=22 admixture run
write.csv(best_qfile %>% 
            dplyr::rename(samples=isotype,cluster=subpop,frac_cluster=fraction) %>%
            dplyr::group_by(samples) %>%
            dplyr::mutate(max_frac=max(frac_cluster)) %>%
            dplyr::select(samples,cluster,frac_cluster,max_frac,lat,long,geo),file = "../../processed_data/genetic_similarity_and_admixutre/best_k_long_admix_pops.csv")

#get list of non-admixed strains and their subpopulations
admix <- best_NADM_assignments %>% dplyr::select(isotype,cluster=subpop)

#write non-admixed isotypes
admix_out <- admix %>% dplyr::left_join(geo_info, by="isotype") %>%
  dplyr::rename(samples=isotype) %>%
  dplyr::select(samples,cluster,lat,long,geo)

write.table(admix_out,file = "../../processed_data/genetic_similarity_and_admixutre/non_admixed_isotypes.txt",sep = "\t",quote = F,row.names = F)

#get subpopulation fractions and geographic information for every isotype strain
pops <- best_qfile %>%
  dplyr::left_join(admix %>% dplyr::rename(admix=cluster),by=c("isotype")) 

#set subpopulation colors, K=22
admix_color <- data.frame(
  letter = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J",
             "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U","V","W","X","Admixed","Cosmopolitan"),
  color = c("#4B0401", "#DA030F", "#FF6C93", "#D45602", "#563901",
            "#FFE5CA", "#FFB914", "#FFDA90", "#77C002", "#05E51B",
            "#02491E", "#05DEA2", "#82FFFA", "#01A4B1", "#5CB7FF",
            "#000F2D", "#0340B9", "#9C87FF", "#CBB1FF", "#F479FF", "#5F0158","#BAB465","#FF007F","#E5A995","gray80","gray30"),
  stringsAsFactors = FALSE
)

#read pairwise similarity estimates
conc <- readr::read_tsv("../../processed_data/genetic_similarity_and_admixutre/gtcheck.tsv")

#construct symmetric pairwise similarity matrix
concordance_matrix <- conc %>%
  dplyr::filter(i %in% isos$isotype & j %in% isos$isotype) %>%
  dplyr::mutate(concordance = (sites - discordance) / sites) %>%
  dplyr::select(i, j, concordance) %>%
  dplyr::bind_rows(
    conc %>%
      dplyr::filter(i %in% isos$isotype & j %in% isos$isotype) %>%
      dplyr::mutate(concordance = (sites - discordance) / sites) %>%
      dplyr::select(i = j, j = i, concordance)
  ) %>%
  tidyr::pivot_wider(names_from = i, values_from = concordance) %>%
  tibble::column_to_rownames("j") %>%
  as.matrix()
all_iso <- sort(unique(lineages$isotype))
concordance_matrix <- concordance_matrix[all_iso, all_iso]
diag(concordance_matrix) <- 1 #pairwise comparisons to self are ommitted in gtcheck, we set to 1

#merge admixture, geographic, and color data into one dataframe
geo <- geo_info %>%
  dplyr::left_join(admix,by="isotype") %>%
  dplyr::mutate(cluster=ifelse(is.na(cluster),"Admixed",cluster)) %>%
  dplyr::left_join(admix_color,by=c("cluster"="letter")) %>%
  dplyr::rename(subpop=cluster)  %>% 
  dplyr::mutate(subpop=ifelse(geo=="Cosmopolitan","Cosmopolitan",subpop)) 

lineage_levels <- c("Tropical","TD1","TH","KD","TD2","TD3","Temperate","ID","NWD","Hubei","Quebec","AD")

#get target admixture run | k

getADM_plots <- function(best_kval,cvmat_long,qlist,admix,admix_color,window_size) {
  plotlist <- list()
  k=1
  for (i in (best_kval-window_size):(best_kval+window_size)) {
    
    target_best_seed <- cvmat_long %>% dplyr::filter(K==i) %>% dplyr::filter(cv_error==min(cv_error))
    target_qfile <- qlist[[i-1]] %>% dplyr::filter(seed==(target_best_seed %>% dplyr::pull(seed)))
    target_pops <- target_qfile %>% dplyr::left_join(admix %>% dplyr::rename(admix=cluster),by=c("isotype"))  
    
    adm_plot_df <- target_pops %>% dplyr::filter(fraction>0.001) %>%
      dplyr::mutate(Lineage=ifelse(Lineage=="Montreal","Quebec",Lineage)) %>%
      dplyr::left_join(admix_color, by = c("admix" = "letter"))
    
    #max number of clusters per sample (for padding)
    Kmax <- adm_plot_df %>%
      dplyr::count(isotype, name = "k") %>%
      dplyr::summarise(Kmax = base::max(k), .groups = "drop") %>%
      dplyr::pull(Kmax)
    
    #build an ordering table: Lineage -> dominant cluster -> sorted fractions
    order_tbl <- adm_plot_df %>%
      dplyr::group_by(isotype) %>%
      dplyr::summarise(
        Lineage = dplyr::first(Lineage),
        dominant_subpop = subpop[base::which.max(fraction)],
        sorted_fracs = base::list(base::sort(fraction, decreasing = TRUE)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        # pad sorted fractions to Kmax with zeros for lexicographic ties
        sorted_fracs = purrr::map(sorted_fracs, ~{
          x <- .x
          base::length(x) <- Kmax
          x[base::is.na(x)] <- 0
          stats::setNames(x, base::paste0("f", base::seq_along(x)))
        }),
        Lineage = base::factor(Lineage, levels = lineage_levels)
      ) %>%
      tidyr::unnest_wider(sorted_fracs) %>%
      dplyr::arrange(
        Lineage,                               # 1) by Lineage
        dominant_subpop,                       # 2) by dominant subpop within Lineage
        dplyr::across(tidyselect::starts_with("f"), dplyr::desc)  # 3) by max, 2nd, 3rd, ...
      )
    
    sample_levels <- order_tbl$isotype
    
    #facet by Lineage for clear separation along the x-axis.
    lineage_facets <- adm_plot_df %>%
      dplyr::mutate(Lineage = base::factor(Lineage, levels = lineage_levels)) %>%
      dplyr::distinct(Lineage) %>%
      dplyr::arrange(Lineage) %>%
      dplyr::pull(Lineage) %>%
      base::as.character()
    
    facet_counts <- adm_plot_df %>%
      dplyr::distinct(Lineage, isotype) %>%
      dplyr::mutate(Lineage = base::factor(Lineage, levels = lineage_levels)) %>%
      dplyr::count(Lineage, name = "n") %>%
      # ensure we have a row for every facet in order
      tidyr::complete(Lineage = base::factor(lineage_facets, levels = lineage_facets), fill = list(n = 0)) %>%
      dplyr::arrange(Lineage) %>%
      dplyr::pull(n)
    
    pleg <- adm_plot_df %>%
      #dplyr::filter(!(Lineage %in% c("TD2","TD3","NWD","Hubei","Quebec","ID"))) %>%
      dplyr::mutate(
        isotype = base::factor(isotype, levels = sample_levels),
        Lineage = base::factor(Lineage, levels = lineage_levels)
      ) %>%
      dplyr::group_by(isotype) %>%
      dplyr::arrange(dplyr::desc(fraction), .by_group = TRUE) %>%
      dplyr::ungroup() %>%
      ggplot(aes(x = isotype, y = fraction, fill = subpop)) +
      geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = stats::setNames(admix_color$color, admix_color$letter)) +
      labs(x = "Sample", y = "Fraction", fill = "ADMIXTURE subpopulation") +
      theme_bw() +
      theme(
        axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank(),
        panel.spacing.x = grid::unit(10, "pt"),   # no horizontal gap
        panel.spacing.y = grid::unit(0, "pt"),
        strip.background = element_blank(),
        strip.placement  = "outside",
        strip.text = element_text(angle = 90, hjust = 0.5,
                                  margin = margin(0, 10, 0, 10, unit = "pt"),lineheight = 20,),
        plot.margin = grid::unit(c(5, 5, 5, 5), "pt"),
        axis.title  = element_blank(),
        axis.ticks  = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.position = "right",        
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.justification.bottom = "right",
        legend.title.position = "top"
      ) +
      coord_cartesian(clip = "off") +
      scale_y_continuous(expand = c(0,0)) +
      ggh4x::facet_grid2(~ Lineage, scales = "free_x", space = "free_x") +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE))
    
    p1 <- adm_plot_df %>%
      dplyr::filter(!(Lineage %in% c("KD","TD2","TD3","NWD","Hubei","Quebec","ID"))) %>%
      dplyr::mutate(
        isotype = base::factor(isotype, levels = sample_levels),
        Lineage = base::factor(Lineage, levels = lineage_levels)
      ) %>%
      dplyr::group_by(isotype) %>%
      dplyr::arrange(dplyr::desc(fraction), .by_group = TRUE) %>%
      dplyr::ungroup() %>%
      ggplot(aes(x = isotype, y = fraction, fill = subpop)) +
      geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = stats::setNames(admix_color$color, admix_color$letter)) +
      labs(x = "Sample", y = "Fraction", fill = "subpop") +
      theme_bw() +
      theme(
        axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank(),
        panel.spacing.x = grid::unit(5, "pt"),   # no horizontal gap
        panel.spacing.y = grid::unit(0, "pt"),
        strip.background = element_blank(),
        strip.placement  = "outside",
        strip.text = element_text(angle = 90, hjust = 0.5,
                                  margin = margin(0, 10, 0, 10, unit = "pt"),lineheight = 20,),
        plot.margin = grid::unit(c(0.5,0.5,0.5,0.5), "pt"),
        axis.title.x  = element_blank(),
        axis.ticks  = element_blank(),
        axis.text = element_blank(),
        legend.position = "none"
      ) +
      coord_cartesian(clip = "off") +
      scale_y_continuous(expand = c(0,0)) +
      ggh4x::facet_grid2(~ Lineage, scales = "free_x", space = "free_x")
    
    p2 <- adm_plot_df %>%
      dplyr::filter((Lineage %in% c("KD","TD2","TD3","NWD","Hubei","Quebec","ID"))) %>%
      dplyr::mutate(
        isotype = base::factor(isotype, levels = sample_levels),
        Lineage = base::factor(Lineage, levels = lineage_levels)
      ) %>%
      dplyr::group_by(isotype) %>%
      dplyr::arrange(dplyr::desc(fraction), .by_group = TRUE) %>%
      dplyr::ungroup() %>%
      ggplot(aes(x = isotype, y = fraction, fill = subpop)) +
      geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = stats::setNames(admix_color$color, admix_color$letter)) +
      labs(x = "Sample", y = "Fraction", fill = "subpop") +
      theme_bw() +
      theme(
        axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank(),
        panel.spacing.x = grid::unit(5, "pt"),   # no horizontal gap
        panel.spacing.y = grid::unit(0, "pt"),
        strip.background = element_blank(),
        strip.placement  = "outside",
        strip.text = element_text(angle = 90, hjust = 0.5,
                                  margin = margin(0, 10, 0, 10, unit = "pt"),lineheight = 30),
        plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "pt"),
        axis.title  = element_blank(),
        axis.ticks  = element_blank(),
        axis.text = element_blank(),
        legend.position = "none"
      ) +
      coord_cartesian(clip = "off") +
      scale_y_continuous(expand = c(0,0)) +
      ggh4x::facet_grid2(~ Lineage, scales = "free_x", space = "free_x")
    
    legend <- cowplot::get_legend(pleg)
    if (k==1) {
      plotlist[[k]] <- cowplot::plot_grid(p1 + ylab(paste0("K=",i)),p2,nrow=1,rel_widths = c(1.5,1),align = "h",axis = "bt")
    } else if (k<(1+2*window_size)) {
      plotlist[[k]] <- cowplot::plot_grid(p1 + theme(strip.text = element_blank()) + ylab(paste0("K=",i)),p2+ theme(strip.text = element_blank()),nrow=1,rel_widths = c(1.5,1),align = "h",axis = "bt")
    } else {
      plotlist[[k]] <- cowplot::plot_grid(p1 + theme(strip.text = element_blank()) + ylab(paste0("K=",i)),p2+ theme(strip.text = element_blank()),nrow=1,rel_widths = c(1.5,1),align = "h",axis = "bt")
      plotlist[[k+1]] <- legend
    }
    k=k+1
  }
  return(plotlist)
}

alladm_plots <- getADM_plots(best_kval,cvmat_long,qlist,admix,admix_color,2)
admplots <- cowplot::plot_grid(alladm_plots[[1]],alladm_plots[[2]],alladm_plots[[3]],alladm_plots[[4]],alladm_plots[[5]],cowplot::plot_grid(NULL,alladm_plots[[6]],nrow=1),nrow=6,ncol=1,align="v",axis="lr",rel_heights = c(1.1,0.7,0.7,0.7,0.7,0.5))

lineage_map <- lineages %>% select(isotype, Lineage)

concord_df <- as.data.frame(as.table(concordance_matrix))
colnames(concord_df) <- c("isotype1", "isotype2", "similarity")

concord_df <- concord_df %>%
  dplyr::left_join(lineage_map, by = c("isotype1" = "isotype")) %>%
  dplyr::rename(Lineage1 = Lineage) %>%
  left_join(lineage_map, by = c("isotype2" = "isotype")) %>%
  dplyr::rename(Lineage2 = Lineage) %>%
  dplyr::filter(isotype1 != isotype2)

# Desired lineage order (alphabetical; change to your desired order if needed)
lvl <- sort(unique(lineages$Lineage))

# Normalize orientation: always store (x,y) so x-index <= y-index in lvl
concord_df_norm <- concord_df %>%
  dplyr::mutate(
    i = match(Lineage1, lvl),
    j = match(Lineage2, lvl),
    Lx = if_else(i <= j, Lineage1, Lineage2),
    Ly = if_else(i <= j, Lineage2, Lineage1)
  ) %>%
  dplyr::select(Lineage1 = Lx, Lineage2 = Ly, similarity)

# Summarize mean similarity for each (Lineage1, Lineage2) in the UPPER triangle
lineage_summary <- concord_df_norm %>%
  dplyr::group_by(Lineage1, Lineage2) %>%
  dplyr::summarise(mean_similarity = mean(similarity, na.rm = TRUE), .groups = "drop")

# Build plotting df with consistent factor levels; reverse Y so upper-tri is above diagonal
keep_levels <- setdiff(lvl, "Quebec")  # drop Quebec if desired
heatmap_long <- lineage_summary %>%
  dplyr::filter(Lineage1 %in% keep_levels, Lineage2 %in% keep_levels) %>%
  dplyr::mutate(
    Lineage1 = factor(Lineage1, levels = keep_levels),
    Lineage2 = factor(Lineage2, levels = rev(keep_levels))
  )

cc_sum_heatmap <- 
  ggplot(heatmap_long, aes(x = Lineage1, y = Lineage2, fill = mean_similarity)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", mean_similarity)), size = 3) +
  scale_fill_gradientn(colours = c("#4575B4", "#87C6C2", "#FFFFE0", "#F4D166", "#D73027"),
                       na.value = "grey90") +
  scale_x_discrete(position = "bottom", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_blank(),
        panel.border = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.95,0.85),
        panel.grid = element_blank()) +
  labs(fill = "Genetic\nsimilarity")
 
ggsave(cc_sum_heatmap,filename = "../../figures/EDF4_concordance_heatmap.png",width = 7,height = 7,units = "in",device = "png",dpi = 600,bg = "white")


annotation_df <-as.data.frame(geo %>% dplyr::select(geo,abslat,subpop,Lineage, Sublineage,lineage_color,sublineage_color))
annotation_df$abslat <- as.numeric(annotation_df$abslat)
rownames(annotation_df) <- geo$isotype

# Assume annotation_df$abslat exists and is numeric
abslat_vec <- annotation_df$abslat
names(abslat_vec) <- rownames(annotation_df)
abslat_vec <- abslat_vec[rownames(concordance_matrix)]  # align order

geo_vec <- annotation_df$geo
names(geo_vec) <- rownames(annotation_df)
geo_vec <- geo_vec[rownames(concordance_matrix)] 

adm_vec <- annotation_df$subpop
names(adm_vec) <- rownames(annotation_df)
adm_vec <- adm_vec[rownames(concordance_matrix)] 

subpop_levels <- setdiff(unique(adm_vec), c("Admixed","Cosmopolitan"))  # all clusters except admixed
adm_vec <- factor(adm_vec, levels = c(sort(subpop_levels),"Admixed", "Cosmopolitan"))

admix_col_vec <- admix_color$color
names(admix_col_vec) <- admix_color$letter

lin_vec <-annotation_df$Lineage
names(lin_vec) <- rownames(annotation_df)
lin_vec <- lin_vec[rownames(concordance_matrix)]

lin_cols <- annotation_df %>% dplyr::select(Lineage,lineage_color) %>% dplyr::distinct(Lineage,.keep_all = T)
lin_col_vec <- lin_cols$lineage_color
names(lin_col_vec) <- lin_cols$Lineage

sublin_rf <- annotation_df %>% dplyr::mutate(Sublineage=ifelse(Sublineage %in% c("TS1","TS2"),Sublineage,NA))  

sublin_cols <- sublin_rf %>% dplyr::select(Sublineage,sublineage_color) %>% dplyr::distinct(Sublineage,.keep_all = T) %>%
  dplyr::filter(Sublineage %in% c("TS1","TS2"))
sublin_col_vec <- sublin_cols$sublineage_color 
names(sublin_col_vec) <- sublin_cols$Sublineage

sublin_vec <- sublin_rf$Sublineage
names(sublin_vec) <- rownames(sublin_rf)
sublin_vec <- sublin_vec[rownames(concordance_matrix)]

lat_col_fun <- circlize::colorRamp2(
  seq(0, 60, length.out = 5),
  c("#D73027", "#EDC948" , "#FFFFE0", "#4575B4", "#000080"),
)

geo_vec_filtered <- factor(geo_vec, levels = c(setdiff(geo_vec, c("Cosmopolitan", "unknown")),c("Cosmopolitan", "unknown")))
lin_vec_filtered <- factor(lin_vec, levels = c("Tropical","TD1","TH","TD2","TD3","Temperate","Quebec","Hubei","ID","NWD","KD","AD"))
# Create row annotation
row_annot <- columnAnnotation(
  `Geo.` = geo_vec_filtered,
  Group = lin_vec_filtered,
  col = list(
             `Geo.` = geo_colors,
             Group= lin_col_vec),
  annotation_legend_param = list(
    Group = list(title= "Relatedness\ngroup", ncol = 2, fontsize = 14, title_gp = grid::gpar(fontsize = 9), labels_gp = grid::gpar(fontsize = 8)),
    `Geo.` = list(title = "Geographic\nregion", ncol = 2,   title_gp = grid::gpar(fontsize = 9), labels_gp = grid::gpar(fontsize = 8))
  ),
  annotation_name_side = "left",
  na_col = "white"
)

bottom_annot <- rowAnnotation(
  `Abs.Lat.` = abslat_vec,
  col = list(`Abs.Lat.` = lat_col_fun),
  annotation_legend_param = list(
    `Abs.Lat.` = list(title = "Absolute\nlatitude", ncol = 1,title_gp = grid::gpar(fontsize = 9), labels_gp = grid::gpar(fontsize = 8))
  ),
  annotation_name_side = "top"
)

# Heatmap color function (5-color scale)
phylo_vals <- as.vector(concordance_matrix)
phylo_vals <- phylo_vals[!is.na(phylo_vals)]  # remove NAs

# Set the breakpoints manually
breaks <- c(
  min(phylo_vals),
  quantile(phylo_vals, 0.25),
  median(phylo_vals),
  quantile(phylo_vals, 0.75),
  1
)
# Define the color mapping with white at the median
phylo_col_fun <- circlize::colorRamp2(
  breaks,
  c("#4575B4", "#87C6C2", "#FFFFE0","#F4D166","#D73027")
)
strain_props <- setNames(isos$count, isos$isotype)
strain_props <- strain_props[colnames(concordance_matrix)]

heatmap_grob <- grid::grid.grabExpr({
  ComplexHeatmap::draw(
    ComplexHeatmap::Heatmap(
      concordance_matrix,
      name = "Genetic\nsimilarity",
      col = phylo_col_fun,
      cluster_rows = T,
      cluster_columns = T,
      clustering_distance_rows = function(m) as.dist(1 - m),
      clustering_distance_columns = function(m) as.dist(1 - m),
      clustering_method_rows = "average",
      clustering_method_columns = "average",
      show_row_names = F,
      show_column_names = F,
      right_annotation = bottom_annot,
      bottom_annotation = row_annot,
      heatmap_legend_param = list(ncol = 2,
                                  title_gp = grid::gpar(fontsize = 9),
                                  labels_gp = grid::gpar(fontsize = 8)),
    ),
    merge_legend = F,
    heatmap_legend_side = "right",
    annotation_legend_side = "bottom"
  )
})

ggmap <- ggplotify::as.ggplot(heatmap_grob)

ggsave(plot = ggmap, filename = "../../figures/Figure2_heatmap_cc_byGeoLat_20251014.png", width = 7, height = 7,bg = "white",device = "png",units = "in",dpi = 600)
