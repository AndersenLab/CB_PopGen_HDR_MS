rm(list = ls())

library(ggplot2)
library(dplyr)
library(maps)

Cb_raw_csv <- readr::read_tsv("../../data/master_sheet_cb/Cb_1972_strains_data.tsv") %>%
  filter(!(isotype %in% c("MY681", "ECA1146", "JU356", "ECA1503"))) %>%
  mutate(
    longitude = ifelse(
      strain %in% c("PB857", "PB858", "PB859"),
      "-84.0554",
      longitude
    ),
    longitude = as.numeric(longitude),
    latitude  = as.numeric(latitude)
  )

NIC174_info <- Cb_raw_csv %>%
  select(strain, isotype, latitude, longitude) %>%
  filter(isotype == "NIC174") %>%
  filter(!is.na(latitude), !is.na(longitude))

world <- map_data("world")

europe_map <- world %>%
  filter(
    region != "Antarctica"
  )

p_map_NIC174_Europe <- ggplot() +
  geom_polygon(
    data = europe_map,
    aes(x = long, y = lat, group = group),
    fill = "gray85",
    color = "gray95",
    linewidth = 0.15
  ) +
  geom_point(
    data = NIC174_info,
    aes(x = longitude, y = latitude),
    shape = 16,
    size = 2,
    color = "#53886C",
    alpha = 0.8
  ) +
  coord_quickmap(
    xlim = c(-28, 55),
    ylim = c(36, 68),
    expand = FALSE
  ) +
  theme_void() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    plot.margin = margin(5, 5, 5, 5)
  )

ggsave(
  "../../figures/FigureS2_NIC174_Europe_map.pdf",
  plot = p_map_NIC174_Europe,
  width = 7.5,
  height = 5,
  units = "in"
)

