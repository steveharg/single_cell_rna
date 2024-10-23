# - function for producing a heatmap annotated with two categories
require(ComplexHeatmap)
require(circlize)

heatmap_two_cat <- function(heat,
                            metadata_column,
                            heat_quantity_name,
                            rows_title,
                            cols_title,
                            row_names_fontsize) {
  # colours <- list("Time" = c("pre"="royalblue", "post"="red2"))
  colours <- list("Time" = c("pre" = "black", "post" = "grey"))
  
  #Set colour palatte
  # myCol <- colorRampPalette(c("dodgerblue", "black", "yellow"))(100)
  myCol <- colorRampPalette(c("dodgerblue", "red2", "yellow"))(100)
  
  heat_range = range(heat)
  myBreaks <- seq(heat_range[1], heat_range[2], length.out = 100)
  # myBreaks <- seq(heat_range[1], 2*heat_range[2]/3, length.out=100)
  
  #Set annotation
  ann <- data.frame(metadata_column)
  colnames(ann) <- "Time"
  colAnn <- HeatmapAnnotation(
    df = ann,
    which = "col",
    col = colours,
    annotation_width = unit(c(2, 2, 2, 2), "cm"),
    gap = unit(0.75, "mm"),
    annotation_legend_param = list(
      Time = list(
        nrow = 3,
        title = "Time",
        title_position = "topcenter",
        legend_direction = "vertical",
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 12, fontface = "bold")
      )
    )
  )
  
  boxplotCol <- HeatmapAnnotation(
    boxplot = anno_boxplot(
      data.matrix(heat),
      border = FALSE,
      gp = gpar(fill = "#CCCCCC"),
      pch = ".",
      size = unit(2, "mm"),
      axis = TRUE,
      axis_param = list(gp = gpar(fontsize =
                                    12), side = "left")
    ),
    annotation_width = unit(c(2.0), "cm"),
    which = "col"
  )
  boxplotRow <- HeatmapAnnotation(
    boxplot = row_anno_boxplot(
      data.matrix(heat),
      border = FALSE,
      gp = gpar(fill =
                  "#CCCCCC"),
      pch = ".",
      size = unit(2, "mm"),
      axis = TRUE,
      axis_param = list(gp = gpar(fontsize =
                                    12), side = "top")
    ),
    annotation_width = unit(c(2.0), "cm"),
    which = "row"
  )
  
  hmap <- Heatmap(
    heat,
    name = heat_quantity_name,
    col = colorRamp2(myBreaks, myCol),
    
    heatmap_legend_param = list(
      color_bar = "continuous",
      legend_direction = "horizontal",
      legend_width = unit(6, "cm"),
      legend_height = unit(4, "cm"),
      title_position = "topcenter",
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 12, fontface = "bold")
    ),
    
    cluster_rows = TRUE,
    show_row_dend = TRUE,
    row_title = rows_title,
    row_title_side = "right",
    row_title_gp = gpar(fontsize = 14, fontface = "bold"),
    row_title_rot = 90,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = row_names_fontsize, fontface =
                          "bold"),
    row_names_side = "left",
    row_dend_width = unit(25, "mm"),
    
    cluster_columns = TRUE,
    show_column_dend = TRUE,
    column_title = cols_title,
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    column_title_rot = 0,
    show_column_names = FALSE,
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_dend_height = unit(25, "mm"),
    
    clustering_distance_columns = function(x)
      as.dist(1 - cor(t(x))),
    clustering_method_columns = "ward.D2",
    clustering_distance_rows = function(x)
      as.dist(1 - cor(t(x))),
    clustering_method_rows = "ward.D2",
    
    top_annotation = colAnn
  )
  
  return(list(hmap, boxplotRow))
  
}