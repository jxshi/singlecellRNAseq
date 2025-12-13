SeuratPlot <- function(object, groupBy, reduction, point_size = 3, palette = "ggsci::category20c_d3") {
  reduction <- rlang::as_name(rlang::enquo(reduction))
  groupBy   <- rlang::as_name(rlang::enquo(groupBy))

  if (!reduction %in% names(object@reductions)) {
    stop(sprintf("Reduction '%s' not found in object@reductions.", reduction))
  }

  if (!groupBy %in% colnames(object@meta.data)) {
    stop(sprintf("Metadata column '%s' not found in object@meta.data.", groupBy))
  }

  embeddings <- object@reductions[[reduction]]@cell.embeddings

  if (ncol(embeddings) < 2) {
    stop(sprintf("Reduction '%s' must contain at least two dimensions for plotting.", reduction))
  }

  plot_data <- embeddings %>%
    as.data.frame() %>%
    dplyr::mutate(Barcode = rownames(.)) %>%
    dplyr::inner_join(object@meta.data %>% tibble::rownames_to_column("Barcode"), by = "Barcode") %>%
    tibble::column_to_rownames("Barcode") %>%
    dplyr::select(1, 2, tidyselect::all_of(groupBy)) %>%
    dplyr::rename(group = tidyselect::all_of(groupBy))

  axis_labels <- colnames(plot_data)[1:2]

  centroids <- stats::aggregate(as.matrix(plot_data[, 1:2]) ~ group, data = plot_data, FUN = mean)

  ggplot(data = plot_data, mapping = aes(x = .data[[axis_labels[1]]], y = .data[[axis_labels[2]]])) +
    geom_point(aes(fill = group), color = "white", shape = 21, stroke = 0.5, size = point_size) +
    geom_point(data = centroids, aes(x = .data[[axis_labels[1]]], y = .data[[axis_labels[2]]], fill = group),
               color = "white", shape = 21, stroke = 0.5, size = point_size) +
    theme_bw() +
    ggrepel::geom_label_repel(
      data    = centroids,
      mapping = aes(label = group, color = group),
      size = 7,
      fill = ggplot2::alpha("white", 0.6),
      fontface = "bold",
      family = "Arial",
      show.legend = FALSE
    ) +
    scale_fill_manual(values = ggplot2::alpha(paletteer::paletteer_d(palette), 0.65)) +
    scale_color_manual(values = ggplot2::alpha(paletteer::paletteer_d(palette), 0.65)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(family = "Arial", size = 12, color = "black"),
      axis.text.y = element_text(family = "Arial", size = 12, color = "black"),
      axis.title.x = element_text(family = "Arial", size = 15, color = "black"),
      axis.title.y = element_text(family = "Arial", size = 15, color = "black"),
      axis.ticks = element_line(color = "black", lineend = "round"),
      legend.position = "bottom",
      legend.text = element_text(family = "Arial", size = 10, color = "black"),
      legend.title = element_text(family = "Arial", size = 13, color = "black")
    )
}
