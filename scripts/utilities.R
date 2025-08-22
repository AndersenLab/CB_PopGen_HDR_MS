
# region colors
geo.colours <- c("Hawaii"="#66C2A5", "Australia"="#FC8D62", "Central America"="#8DA0CB",
                 "South America"="#E78AC3", "Africa"="#A6D854", "Caribbean"="#FFD92F",
                 "Taiwan" = "#E5C494", "North America" = "#A65628", "Europe" = "#E41A1C",
                 "Asia" = "#399EB8", "Pacific" = "#611EA1",
                 "New Zealand" = "green","Atlantic" = "purple",
                 "Oceania" ="#DB7779",
                 "Micronesia" = "#E7298A",
                 "Indonesia" = "#7570B3",
                 "Malay Archipelago" = "#4110B3",
                 "C. America"="#8DA0CB",
                 "S. America"="#E78AC3",
                 "N. America" = "#A65628",
                 # "Globalized" = "gray30",
                 "cosmopolitan" = "gray30",
                 "unknown" = 'grey'
                 )


species_color <- c("C. elegans"   = "#DB6333",
                   "C. tropicalis" = "#0719BC",
                   "C. briggsae"   = "#53886C")


genome_domain_colors <- c("Tip" = "#5E3C99", 
                          "Center" = "#FDB863", 
                          "Arm" = "#4393C3")


library(pals)
# Cluster color palettes 
group_col <- kelly()[c(3, 4, 5, 6, 7, 8, 10, 9)]
names(group_col) <- c("1", "2", "3", "4", "5", "6", "7", "outlier")
group_col8 <- kelly()[c(3, 4, 5, 6, 7, 8, 10, 11, 9)]
names(group_col8) <- c("1", "2", "3", "4", "5", "6", "7", "8", "outlier")
group_col9 <- kelly()[c(3, 4, 5, 6, 7, 8, 10, 11, 12, 9)]
names(group_col9) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "outlier")


# Hawaii island color palettes
isl_cols <- alphabet()[c(1, 2, 3, 4, 6, 5)]
names(isl_cols) <- c("Kauai", "Oahu", "Molokai", "Maui", "Big Island", "unknown")







# Function to grab ggrepel boxes
findboxes <- function(
    df, xcol, ycol,
    box_padding_x, box_padding_y,
    point_padding_x, point_padding_y,
    xlim, ylim,
    force = 1e-7, maxiter = 20000
) {
  
  # x and y posiitons as a dataframe
  posdf <- df[c(xcol, ycol)]
  
  # returnd a df where columns are points
  boxdf <- apply(posdf, 1, function(row) {
    xval <- row[xcol]
    yval <- row[ycol]
    return(c(
      xval - box_padding_x / 2,
      yval - box_padding_y / 2,
      xval + box_padding_x / 2,
      yval + box_padding_y / 2
    ))
  })
  # columns are x1,y1,x2,y2
  boxmatrix <- as.matrix(t(boxdf))
  
  moved <- ggrepel:::repel_boxes(
    data_points = as.matrix(posdf),
    point_padding_x = point_padding_x,
    point_padding_y = point_padding_y,
    boxes = boxmatrix,
    xlim = xlim,
    ylim = ylim,
    hjust = 0.5,
    vjust = 0.5,
    force = force,
    maxiter = maxiter
  )
  
  finaldf <- cbind(df$isotype, posdf, moved)
  names(finaldf) <- c("isotype", "x1", "y1", "x2", "y2")
  return(finaldf)
}



# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}









