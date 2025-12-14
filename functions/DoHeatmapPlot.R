DoHeatmapPlot <- function(object, groupBy, features) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required for DoHeatmapPlot().")
  }
  if (!requireNamespace("BuenColors", quietly = TRUE)) {
    stop("Package 'BuenColors' is required for DoHeatmapPlot().")
  }
  if (missing(groupBy) || !is.character(groupBy) || length(groupBy) != 1) {
    stop("`groupBy` must be a single metadata column name.")
  }
  if (missing(features) || length(features) == 0) {
    stop("`features` must contain at least one gene name.")
  }

  meta_cols <- colnames(object@meta.data)
  if (!groupBy %in% meta_cols) {
    stop("`groupBy` must be a column in object@meta.data.")
  }

  available_features <- rownames(object)
  if (any(!features %in% available_features)) {
    missing_feats <- setdiff(features, available_features)
    stop("The following features are missing from the object: ", paste(missing_feats, collapse = ", "))
  }

  plot_data <- SeuratObject::FetchData(
    object = object,
    vars = c(features, groupBy),
    slot = "counts"
  ) %>%
    dplyr::mutate(across(.cols = where(is.numeric), .fns = ~ log2(.x + 1))) %>%
    dplyr::rename(group = as.name(groupBy)) %>%
    dplyr::arrange(group)

  clusterInfo <- plot_data$group
  if (length(unique(clusterInfo)) == 0) {
    stop("No groups found to split columns by.")
  }

  plot_matrix <- plot_data %>%
    dplyr::select(-group) %>%
    t()

  col <- BuenColors::jdb_color_maps[seq_along(unique(clusterInfo))]
  names(col) <- as.character(unique(clusterInfo))

  top_anno <- ComplexHeatmap::HeatmapAnnotation(
    cluster = ComplexHeatmap::anno_block(
      gp = grid::gpar(fill = col),
      labels = as.character(unique(clusterInfo)),
      labels_gp = grid::gpar(cex = 0.5, col = "white", family = "Arial")
    )
  )

  plot_matrix <- as.matrix(as.data.frame(plot_matrix))
  gene_pos <- which(rownames(plot_matrix) %in% features)

  row_anno <- ComplexHeatmap::rowAnnotation(
    mark_gene = ComplexHeatmap::anno_mark(at = gene_pos, labels = features)
  )

  ComplexHeatmap::Heatmap(
    matrix = plot_matrix,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    show_heatmap_legend = TRUE,
    column_split = clusterInfo,
    top_annotation = top_anno,
    column_title = NULL,
    right_annotation = row_anno,
    use_raster = FALSE,
    heatmap_legend_param = list(title = "log2(count+1)", title_position = "leftcenter-rot")
  )
}
