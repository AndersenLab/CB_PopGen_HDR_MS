library(ape)
library(ggplot2)
library(ggtree)
library(cowplot)
library(dplyr)
library(tidyr)

mastsumm <- readr::read_tsv("../../processed_data/tree_topology/mast/mast_summary.tsv") #%>%
matrixsumm <- readr::read_tsv("../../processed_data/tree_topology/phylo_py_pruner/supermatrix_stats_summary.tsv") %>%
  dplyr::mutate(`Gene set`=ifelse(`Gene set`=="ALL","mast_all","mast_hdr"))
quartetsumm <- base::list.files(
  "../../processed_data/tree_topology/mast/quartet_topologies/summ_quartet_topos",
  pattern = "_quartets_summary\\.tsv$",
  full.names = TRUE
) %>%
  rlang::set_names() %>%
  purrr::map_dfr(
    function(.x) {
      readr::read_tsv(.x, show_col_types = FALSE) %>%
        dplyr::mutate(
          CBY = stringr::str_extract(base::basename(.x), "^[^_]+")
        )
    }
  ) %>%
  dplyr::mutate(topology=gsub("CB2","CBY",topology),
                topology=gsub("CB1","CBX",topology)) %>%
  dplyr::rename(tree=quartet) %>%
  dplyr::select(-samples) %>%
  dplyr::group_by(quartet_id, CBY) %>%
  dplyr::mutate(
    topo_num = paste0("T", dplyr::row_number())
  ) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(
    id_cols = c(quartet_id, CBY),
    names_from = topo_num,
    values_from = c(topology, tree),
    names_glue = "{.value}_{topo_num}"
  ) %>%
  dplyr::rename(
    T1 = topology_T1,
    T2 = topology_T2,
    T3 = topology_T3,
    Tree1 = tree_T1,
    Tree2 = tree_T2,
    Tree3 = tree_T3
  )

comb_summ <- mastsumm %>% 
  dplyr::left_join(matrixsumm,by=c("quartet"="Quartet ID","set"="Gene set","tag"="CB Background")) %>%
  dplyr::left_join(quartetsumm,by=c("quartet"="quartet_id","tag"="CBY"),relationship = "many-to-many")


bests1_byBIC <- comb_summ %>% 
  dplyr::group_by(set,tag,quartet) %>%
  dplyr::filter(BIC==min(BIC)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(set) %>%
  dplyr::summarise(
    n_best_sub1 = sum(submodel == "submodel1"),
    .groups = "drop"
  )

write.table(comb_summ %>% dplyr::select(-occ,-status),"../../supplementary_data/ST12_mast_summary.tsv",quote = F,row.names = F,sep = "\t")

#ggplot(comb_summ %>% dplyr::filter(set=="mast_hdr")) + geom_point(aes(x=TW1,y=`Number of alignments`))

mastsumm_long <- mastsumm %>% 
  dplyr::filter(submodel=="submodel1") %>%
  dplyr::filter(set %in% c("mast_hdr", "mast_all")) %>%
  tidyr::pivot_longer(
    cols = c(TW1, TW2, TW3),
    names_to = "tree",
    values_to = "weight"
  ) %>%
  dplyr::mutate(
    set = dplyr::recode(
      set,
      mast_hdr = "HDR genes",
      mast_all = "All genes"
    ),
    occ = dplyr::recode(
      occ,
      occ80 = "80% Occupancy"
    ),
    tree = factor(tree, levels = c("TW1", "TW2", "TW3"))
  ) %>%
  dplyr::mutate(
    tag = factor(
      tag,
      levels = c("JU3237", "JU3200", "QG4232")
    )
  )

f5b <- ggplot2::ggplot(
  mastsumm_long,
  ggplot2::aes(x = tree, y = weight)) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  ggplot2::geom_point(
    ggplot2::aes(color = quartet),
    position = ggplot2::position_jitter(width = 0.15, height = 0),
    size = 1,
    alpha = 0.8
  ) +
  ggplot2::facet_grid(set ~ tag) +
  ggplot2::scale_x_discrete(
    labels = c(
      TW1 = expression(T[1]),
      TW2 = expression(T[2]),
      TW3 = expression(T[3])
    )
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  ggplot2::labs(
    x = "Topology",
    y = "Tree weight",
    color = "Quartet"
  )

mastsumm_long6 <- mastsumm %>% 
  dplyr::filter(submodel=="submodel6") %>%
  dplyr::filter(set %in% c("mast_hdr", "mast_all")) %>%
  tidyr::pivot_longer(
    cols = c(TW1, TW2, TW3),
    names_to = "tree",
    values_to = "weight"
  ) %>%
  dplyr::mutate(
    set = dplyr::recode(
      set,
      mast_hdr = "HDR genes",
      mast_all = "All genes"
    ),
    occ = dplyr::recode(
      occ,
      occ80 = "80% Occupancy"
    ),
    tree = factor(tree, levels = c("TW1", "TW2", "TW3"))
  ) %>%
  dplyr::mutate(
    tag = factor(
      tag,
      levels = c("JU3237", "JU3200", "QG4232")
    )
  )

sf21 <- ggplot2::ggplot(
  mastsumm_long6,
  ggplot2::aes(x = tree, y = weight)) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  ggplot2::geom_point(
    ggplot2::aes(color = quartet),
    position = ggplot2::position_jitter(width = 0.15, height = 0),
    size = 1,
    alpha = 0.8
  ) +
  ggplot2::facet_grid(set ~ tag) +
  ggplot2::scale_x_discrete(
    labels = c(
      TW1 = expression(T[1]),
      TW2 = expression(T[2]),
      TW3 = expression(T[3])
    )
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  ggplot2::labs(
    x = "Topology",
    y = "Tree weight",
    color = "Quartet"
  )


minor_topos <- mastsumm %>%
  dplyr::filter(submodel=="submodel1") %>%
  dplyr::filter(set %in% c("mast_all", "mast_hdr")) %>%
  dplyr::mutate(`TW2+TW3`=TW2+TW3) %>%
  tidyr::pivot_longer(
    cols = c(TW2, TW3,`TW2+TW3`),
    names_to = "tree",
    values_to = "weight"
  ) %>%
  dplyr::mutate(
    set = dplyr::recode(
      set,
      mast_hdr = "HDR genes",
      mast_all = "All genes"
    ),
    occ = dplyr::recode(
      occ,
      occ80 = "80% Occupancy"
    )
  ) %>%
  dplyr::mutate(
    tag = factor(
      tag,
      levels = c("JU3237", "JU3200", "QG4232")
    ),
    tree = factor(
      tree,
      levels = c("TW2", "TW3", "TW2+TW3"),
      labels = c("T[2]", "T[3]", "T[2]+T[3]")
    )
  )

sf22 <- ggplot2::ggplot(
  minor_topos,
  ggplot2::aes(x = set, y = weight)
) +
  ggplot2::geom_line(
    ggplot2::aes(group = quartet, color = quartet),
    alpha = 0.6
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = quartet),
    size = 2,
    alpha = 0.9
  ) +
  ggplot2::facet_grid(tree ~ tag,labeller = ggplot2::label_parsed) +
  ggplot2::scale_y_continuous(limits = c(0, NA)) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "bottom"
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Tree weight",
    color = "Quartet"
  )

minor_topos_tree_x <- mastsumm %>%
  dplyr::filter(submodel=="submodel1") %>%
  dplyr::filter(set %in% c("mast_all", "mast_hdr")) %>%
  tidyr::pivot_longer(
    cols = c(TW2, TW3),
    names_to = "tree",
    values_to = "weight"
  ) %>%
  dplyr::mutate(
    set = dplyr::recode(
      set,
      mast_hdr = "HDR genes",
      mast_all = "All genes"
    ),
    occ = dplyr::recode(
      occ,
      occ80 = "80% Occupancy"
    )
  ) %>%
  dplyr::mutate(
    tag = factor(
      tag,
      levels = c("JU3237", "JU3200", "QG4232")
    ),
    tree = factor(tree, levels = c("TW2", "TW3"))
  )

sf23 <- ggplot2::ggplot(
  minor_topos_tree_x,
  ggplot2::aes(x = tree, y = weight)
) +
  ggplot2::geom_line(
    ggplot2::aes(group = quartet, color = quartet),
    alpha = 0.6
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = quartet),
    size = 2,
    alpha = 0.9
  ) +
  ggplot2::facet_grid(set ~ tag) +
  ggplot2::scale_x_discrete(
    labels = c(
      TW2 = expression(T[2]),
      TW3 = expression(T[3])
    )
  ) +
  ggplot2::scale_y_continuous(limits = c(0, NA)) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "bottom"
  ) +
  ggplot2::labs(
    x = "Topology",
    y = "Tree weight",
    color = "Quartet"
  )

# trees <- read.tree(text = c(
#   "(Macaca_fascicularis,(Macaca_mulatta,Colobus_angolensis_palliatus),Macaca_nemestrina);",
#   "(Macaca_fascicularis,(Macaca_mulatta,Macaca_nemestrina),Colobus_angolensis_palliatus);",
#   "(Macaca_fascicularis,Macaca_mulatta,(Macaca_nemestrina,Colobus_angolensis_palliatus));"
# ))

trees <- read.tree(text = c(
  "(C. remanei,(C. briggsae Y,C. briggsae X),C. nigoni);",
  "(C. remanei,(C. nigoni,C. briggsae X),C. briggsae Y);",
  "(C. remanei,(C. nigoni,C. briggsae Y),C. briggsae X);"
))

trees <- lapply(
  trees,
  ape::root,
  outgroup = "C.remanei",
  resolve.root = TRUE
)


plots <- lapply(
  seq_along(trees),
  function(i) {
    
    lab_map <- c(
      "C.remanei"    = "~italic(C.~remanei)",
      "C.briggsaeX"  = "~italic(C.~briggsae)~X",
      "C.briggsaeY"  = "~italic(C.~briggsae)~Y",
      "C.nigoni"     = "~italic(C.~nigoni)"
    )
    
    p <- ggtree::ggtree(trees[[i]])
    
    p$data$species <- dplyr::case_when(
      p$data$label == "C.remanei" ~ "remanei",
      grepl("^C\\.briggsae", p$data$label) ~ "briggsae",
      p$data$label == "C.nigoni" ~ "nigoni",
      TRUE ~ NA_character_
    )
    
    p +
      ggtree::geom_tiplab(
        ggplot2::aes(
          label = lab_map[label],
          color = species
        ),
        parse = TRUE
      ) +
      ggplot2::scale_color_manual(
        values = c(
          remanei  = "forestgreen",
          briggsae = "red",
          nigoni   = "blue"
        )
      ) +
      ggplot2::xlim(0, 7) +
      ggplot2::theme(
        legend.position = "none"
      )
  }
)

f5a <- cowplot::plot_grid(plotlist = plots,nrow=1,labels = c("T₁", "T₂", "T₃"),label_x = 0.5,label_fontface = "plain",label_size = 10)

f5 <- cowplot::plot_grid(f5a,f5b,nrow=2,labels=c("a","b"), rel_heights = c(0.5,1))

ggsave(f5,filename = "../../figures/Figure5_mast_results.png",width = 7,height = 5.5,device = "png",units = "in",bg ="white",dpi = 900)
ggsave(sf21,filename = "../../figures/SF21_mast_s6.png",width = 7,height = 5,device = "png",units = "in",bg ="white",dpi = 900)
ggsave(sf22,filename = "../../figures/SF22_mast_asym.png",width = 7,height = 7,device = "png",units = "in",bg ="white",dpi = 900)
ggsave(sf23,filename = "../../figures/SF23_mast_minor_topo.png",width = 7,height = 5.5,device = "png",units = "in",bg ="white",dpi = 900)


