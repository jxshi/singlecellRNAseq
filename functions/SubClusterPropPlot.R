SubPropPlot <- function(object, groupBy) {
    if (missing(groupBy)) {
        stop("`groupBy` must be provided and refer to a metadata column.")
    }
    required_cols <- c(as.character(substitute(groupBy)), "subcelltype")
    missing_cols <- setdiff(required_cols, colnames(object@meta.data))
    if (length(missing_cols) > 0) {
        stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
    }

    plot_data <- object@meta.data %>%
      dplyr::select({{groupBy}}, subcelltype) %>%
      dplyr::rename(group = as.name(groupBy), seurat_clusters = .data$subcelltype) %>%
      tidyr::drop_na()

    figure <- ggstatsplot::ggbarstats(
      data = plot_data,
      x = seurat_clusters,
      y = group,
      package = "ggsci",
      palette = "category20c_d3",
      results.subtitle = FALSE,
      bf.message = FALSE,
      proportion.test = FALSE,
      label.args = list(
        size = 2,
        fill = "white",
        alpha = 0.85,
        family = "Arial",
        fontface = "bold"
      ),
      perc.k = 2,
      title = "",
      xlab = "",
      legend.title = "Seurat Cluster",
      ggtheme = ggpubr::theme_pubclean()
    ) +
      theme(
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black", lineend = "round"),
        legend.position = "right",
        axis.text.x = element_text(size = 12, color = "black", family = "Arial"),
        axis.text.y = element_text(size = 12, color = "black", family = "Arial"),
        legend.text = element_text(family = "Arial", size = 10, color = "black"),
        legend.title = element_text(family = "Arial", size = 13, color = "black")
      )

    gginnards::delete_layers(x = figure, match_type = "GeomText")
}
