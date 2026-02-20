library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

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

df_colors <- data.frame(unname(geo_colors),names(geo_colors)) %>% dplyr::rename(color=`unname.geo_colors.`, geo=`names.geo_colors.`)
#tree <- ape::read.tree(file="eigenstrat_LD0.7_input.min4.phy.treefile")

#tree <- root(tree, outgroup = "ECA2666", resolve.root = TRUE)
isos <- readr::read_tsv(file="../../processed_data/genetic_similarity_and_admixutre/isotype_groups.tsv") %>%
  dplyr::group_by(isotype) %>%
  dplyr::summarise(count=n())

admix <- readr::read_tsv(file="../../processed_data/genetic_similarity_and_admixutre/non_admixed_isotypes.txt") %>%
  dplyr::select(samples,cluster)

pops <- readr::read_csv("../../processed_data/genetic_similarity_and_admixutre/best_k_long_admix_pops.csv") %>%
  dplyr::left_join(admix %>% dplyr::rename(admix=cluster),by=c("samples"))

#set subpopulation colors
admix_color <- data.frame(
  letter = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J",
             "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U","V","W","X","Admixed","Cosmopolitan"),
  color = c("#4B0401", "#DA030F", "#FF6C93", "#D45602", "#563901",
            "#FFE5CA", "#FFB914", "#FFDA90", "#77C002", "#05E51B",
            "#02491E", "#05DEA2", "#82FFFA", "#01A4B1", "#5CB7FF",
            "#000F2D", "#0340B9", "#9C87FF", "#CBB1FF", "#F479FF", "#5F0158","#BAB465","#FF007F","#E5A995","gray80","gray30"),
  stringsAsFactors = FALSE
)


# tree_dend <- ReadDendrogram(file="/vast/eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/test_GTR_LD09/phy_file_LD_0.9.phy.contree.rooted")
# tree_nwk <- ape::read.tree(file="/vast/eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/test_GTR_LD09/phy_file_LD_0.9.phy.contree.rooted")

geo <- readr::read_csv(file="../../processed_data/Geo_info/Cb_indep_isotype_info_geo.csv") %>%
  dplyr::left_join(df_colors,by=c("geo")) %>%
  dplyr::mutate(abslat=abs(lat)) %>%
  dplyr::left_join(admix,by=c("isotype"="samples")) %>%
  dplyr::mutate(cluster=ifelse(is.na(cluster),"Admixed",cluster)) %>%
  dplyr::left_join(admix_color,by=c("cluster"="letter")) %>%
  dplyr::rename(subpop=cluster)  %>% 
  dplyr::mutate(subpop=ifelse(geo=="Cosmopolitan","Cosmopolitan",subpop))

#a genetic similarity dendrogram was uploaded to itol.embl.de
#we collapsed 12 major "clades" that define the relatedness groups and extracted leaf ID (isotype strain name) lists using iToL built-in options
#we used chatGPT to convert the leaf ID lists into R vectors
#the R vectors were pasted below
NWD <- c("JU2767","ED3102","ECA276","ECA163")

indonesia <- c("HPT18", "HPT24", "HPT11")

hubei <- c("VX34", "JU3326")

quebec <- c("QR25")

kerala <- c("JU1341", "JU1348", "VSL2207", "VSL2209",
            "JU3206", "JU3201", "JU3200", "JU2801",
            "JU3203", "JU3202")

twd3 <- c("NIC893", "BRC20503", "BRC20502")

twd2 <- c("NIC819", "BRC20541", "NIC1607", "BRC20115", "BRC20519", "BRC20520", "BRC20492")

twd1 <- c("BRC20324", "TWN1824", "BRC20232", "NIC827", "BRC20229", "NIC19", "NIC1540", "BRC20343",
          "BRC20524", "BRC20522", "NIC20", "NIC1593", "NIC21", "NIC1589", "TWN1899", "NIC22",
          "BRC20075", "BRC20530", "NIC1665", "NIC1664", "NIC1591", "BRC20347", "BRC20531", "NIC825",
          "NIC816", "NIC828", "BRC20239")

australia <- c("QG2908", "QG2893", "QG4232", "QG4208", "QG4093", "QG2919",
               "ECA2670", "ECA2621", "ECA2639", "QG4031", "QG4026", "ECA2617",
               "QG4028", "ECA2618", "QG2891", "QG4063", "QG4132", "QG4096",
               "QG4131", "QG4097", "QG4225", "QG4214", "QG4104", "QG4064",
               "QG4114", "QG4068", "QG4067", "ECA2666", "ECA2647", "QG4095",
               "QG4233", "QG4094")

central_sub  <- c("NIC96", "NIC1491", "NIC1442", "NIC1135", "NIC1190", "NIC1138", "NIC1124", "NIC1141", "NIC1183", "JU2518",
                  "NIC60", "NIC331", "NIC1421", "NIC330", "NIC333", "NIC400", "NIC332", "NIC327", "NIC402", "NIC329",
                  "NIC325", "NIC309", "NIC1304", "NIC1292", "NIC1302", "NIC1267", "NIC1283", "NIC1407", "NIC1423", "NIC1288",
                  "QG3931", "QG3924", "NIC392", "JU1399", "NIC1103", "NIC1050", "NIC1052", "JU2597", "NIC1054", "NIC1059",
                  "NIC1055", "NIC1057", "NIC1060", "NIC1058", "QG761", "QG860", "QG881", "QG700", "QG795", "QG3056",
                  "QG2966", "QG3004", "QG3032", "QG1155", "QG887", "QG997", "QG865", "QG926", "QG856", "QG2663",
                  "QG866", "QG808", "QG855", "QG791", "QG991", "QG874", "QG825", "QG952", "QG986", "QG948",
                  "QG980", "QG727", "QG1141", "QG947", "QG879", "QG833", "QG781", "QG807", "QG828", "QG788",
                  "QG827", "QG878", "QG806", "QG739", "QG2683", "QG995", "QG802", "QG868", "QG733", "QG2964",
                  "QG763", "QG3058", "QG880", "QG736", "QG776", "QG3019", "QG1005", "QG3867", "QG3801", "QG2736",
                  "QG2659", "QG1111", "QG784", "QG2641", "QG737", "QG2972", "QG978", "QG864", "QG809", "QG3823",
                  "QG790", "QG3006", "QG858", "QG2665", "QG3057", "QG875", "QG804", "QG786", "QG994", "QG888",
                  "QG946", "QG945", "QG797", "QG933", "QG3025", "QG3034", "QG3030", "QG2996", "QG2648", "QG998",
                  "QG1007", "QG3041", "QG1126", "QG937", "QG2969")

temperate <- c("PB800","QG3661", "GUN124", "NIC174", "SOW22", "JU3272", "JU2536", "NIC1632", "SOW18", "NIC899",
               "JU2872", "JU3416", "JU1257", "HK104", "BW287", "SOW21", "JU1564", "QG588", "BRC20388",
               "NIC1635", "NIC1633", "PB857", "ECA569", "QX1547", "ECA278", "QG581", "QG547", "EG4360","JU2457")

trop_islandic_sub <- c("JU2162", "JU2160", "JU3237", "JU2470", "JU2057", "VSL2217", "JU2507", "BRC20257", "QG5589",
                       "QG4647", "NIC791", "BRC20119", "QG4805", "QG4701", "QG5603", "QG5474", "ZZY1123", "ZZY0992",
                       "HPT62", "QG2884", "QG2906", "QG2895", "QG2902", "QG2903", "QG2892", "QG2925", "QG2923", "QG2880")

ht_lineage <- c("QG133", "EG6268", "EG6265", "UH1", "ECA1383", "QX1798", "ECA3502", "ECA1170","ECA1443",
                "BRC20333", "NIC1590", "NIC779", "NIC863", "NIC865", "BRC20557", "BRC20244", "JU2216", "BRC20367",
                "BRC20341", "NIC505", "BRC20529", "NIC1666", "NIC1663", "NIC1661", "BRC20230", "NIC814", "BRC20088",
                "NIC1681", "TWN3077", "BRC20081", "BRC20532", "NIC1205", "BRC20463", "BRC20537", "NIC815", "NIC1220",
                "NIC897", "NIC805", "NIC824", "NIC1595", "NIC1660", "NIC889", "NIC864", "BRC20543", "NIC1662", 
                "ECA2385", "BRC20246", "BRC20238", "NIC1693", "NIC1631", "BRC20390", "NIC1692", "BRC20348", "ECA287",
                "TWN1772", "NIC504", "BRC20339",  "BRC20320", "ECA2299", "BRC20290", "BRC20289", "BRC20272",
                "BRC20270", "BRC20545", "NIC822", "BRC20346", "NIC1645", "BRC20499", "NIC1602", "NIC1598",
                "NIC1630", "NIC1562", "BRC20533", "BRC20536", "BRC20534", "BRC20544", "BRC20535", "BRC20501",
                "BRC20366", "TWN2919", "NIC1667", "BRC20334", "BRC20311", "NIC1539", "NIC519", "NIC518",
                "BRC20331", "TWN1907", "BRC20312",  "BRC20315", "BRC20370", "BRC20309", "BRC20308","ECA1559","ECA1605") 

#we detected abnormally increased heterozygosity in these strains (likely an accidental mixture of two strains during library prep)
#they are flagged here, but later excluded from any analysis
flagged_mixtures  <- c("MY681",
                       "JU356",
                       "ECA1146",
                       "ECA1503")

all_divergent <- c(australia,twd1,twd2,twd3,kerala,quebec,hubei,indonesia,NWD)

vector_list <- list(
  NWD = NWD,
  ID = indonesia,
  Hubei = hubei,
  Quebec = quebec,
  TD1 = twd1,
  TD2 = twd2,
  TD3 = twd3,
  KD = kerala,
  AD = australia,
  `Tropical-Central_sub` = central_sub,
  Temperate = temperate,
  #`Temperate-Adjacent_sub` = temperate_adjacent_sub,
  `Tropical-Islandic_sub` = trop_islandic_sub,
  TH = ht_lineage,
  FM = flagged_mixtures
)

collapsed_df <- do.call(rbind, lapply(names(vector_list), function(name) {
  data.frame(Strain = vector_list[[name]], Group = name, stringsAsFactors = FALSE)
}))

tropical_tropical <- data.frame(Strain=isos$isotype[!(isos$isotype %in% collapsed_df$Strain)],Group="Tropical_sub")
all_iso_byLinage <- rbind(collapsed_df,tropical_tropical) %>%
  dplyr::mutate(Lineage=ifelse(grepl("sub",Group),"Tropical",Group)) %>%
  tidyr::separate(Group,into=c("Sublineage","cat"),sep="_sub") %>%
  dplyr::select(-cat) %>%
  dplyr::left_join(geo,by=c("Strain"="isotype")) %>%
  dplyr::rename(isotype=Strain) %>%
  dplyr::rename(geo_color=color.x,adm_color=color.y) %>%
  dplyr::mutate(adm_color=ifelse(geo=="Cosmopolitan","white",adm_color)) %>%
  dplyr::mutate(Sublineage=ifelse(Sublineage=="Tropical" & Lineage=="Tropical","TS2",Sublineage)) %>%
  #dplyr::mutate(Sublineage=ifelse(Sublineage=="Temperate-Adjacent" & Lineage=="Tropical","TS1",Sublineage)) %>%
  dplyr::mutate(Sublineage=ifelse(Sublineage=="Tropical-Islandic" & Lineage=="Tropical","TS1",Sublineage)) %>%
  dplyr::mutate(Sublineage=ifelse(Sublineage=="Tropical-Central" & Lineage=="Tropical","TS2",Sublineage))


lin_colors <- c(TS2="#ff77ab", TC="#ff0000", Tropical="#ff0000",TS3="#ff0000",NWD="#00ff00",#WD="#d0fe02", 
                Quebec="#ff037e", Hubei="#16537e", ID="#6a329f", TD1="#f3c588", TD2="#ffdd02",
                TD3="#00f4c2", KD = "#8b3700", Temperate="#0000ff", TS1="#9fc5e8", AD="#ff8200", TH="#aeb400",FM="#000000")

lin_colors_df <- data.frame(
  sublineage = names(lin_colors),
  color = as.character(lin_colors),
  stringsAsFactors = FALSE
)

all_iso_byLineage_wCol <- all_iso_byLinage %>% 
  dplyr::left_join(lin_colors_df, by=c("Lineage"="sublineage")) %>%
  dplyr::rename(lineage_color=color) %>%
  dplyr::left_join(lin_colors_df, by=c("Sublineage"="sublineage")) %>%
  dplyr::rename(sublineage_color=color) 


### test colors
lineage_legend_df <- all_iso_byLineage_wCol %>%
  dplyr::select(Lineage, lineage_color) %>%
  dplyr::distinct() %>%
  dplyr::mutate(type = "Lineage") %>%
  dplyr::rename(group = Lineage, color = lineage_color)

# Create unique legend entries for sublineage
sublineage_legend_df <- all_iso_byLineage_wCol %>%
  dplyr::select(Sublineage, sublineage_color) %>%
  dplyr::distinct() %>%
  dplyr::mutate(type = "Sublineage") %>%
  dplyr::rename(group = Sublineage, color = sublineage_color)

# Combine both for plotting
legend_df <- dplyr::bind_rows(lineage_legend_df, sublineage_legend_df)

# Dummy plot for legend color preview
# ggplot(legend_df, ggplot2::aes(x = type, y = group, fill = color)) +
#   geom_tile() +
#   scale_fill_identity() +
#   theme_minimal() +
#   labs(title = "Test Legend for Lineage and Sublineage Colors") +
#   theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

write.table(all_iso_byLineage_wCol, file="../../processed_data/genetic_similarity_and_admixutre/isotype_byRG_GeoLocAdmCol_20250909.tsv",quote = F,sep = "\t",row.names = F)
