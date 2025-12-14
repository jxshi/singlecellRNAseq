PropPlot <- function(object, groupBy, colors = NULL, legend_title = "Seurat Cluster") {
  stopifnot("meta.data slot is empty" = !is.null(object@meta.data))

  if (!"celltype" %in% colnames(object@meta.data)) {
    stop("`celltype` column is required in `object@meta.data` for PropPlot().")
  }

  if (rlang::as_name(rlang::enquo(groupBy)) %in% c("celltype", "")) {
    stop("`groupBy` must refer to a metadata column other than `celltype`.")
  }

  if (!rlang::as_name(rlang::enquo(groupBy)) %in% colnames(object@meta.data)) {
    stop(sprintf("`%s` was not found in `object@meta.data`.", rlang::as_name(rlang::enquo(groupBy))))
  }

  plot_data <- object@meta.data %>%
    dplyr::select({{ groupBy }}, celltype) %>%
    dplyr::rename(group = {{ groupBy }})

  base_theme <- theme(
    axis.ticks.x    = element_blank(),
    axis.ticks.y    = element_line(color = "black", lineend = "round"),
    legend.position = "right",
    axis.text.x     = element_text(size = 15, color = "black", family = "Arial", angle = 45, hjust = 1, vjust = 1),
    axis.text.y     = element_text(size = 15, color = "black", family = "Arial"),
    legend.text     = element_text(family = "Arial", size = 10, color = "black"),
    legend.title    = element_text(family = "Arial", size = 13, color = "black")
  )

  figure <- ggstatsplot::ggbarstats(
    data             = plot_data,
    x                = celltype,
    y                = group,
    package          = "ggsci",
    palette          = "category20c_d3",
    results.subtitle = FALSE,
    bf.message       = FALSE,
    proportion.test  = FALSE,
    label.args       = list(
      size     = 2,
      fill     = "white",
      alpha    = 0.85,
      family   = "Arial",
      fontface = "bold"
    ),
    perc.k        = 2,
    title         = "",
    xlab          = "",
    legend.title  = legend_title,
    ggtheme       = ggpubr::theme_pubclean()
  ) + base_theme

  if (!is.null(colors)) {
    if (length(colors) < length(unique(plot_data$celltype))) {
      stop("`colors` must contain at least as many values as there are unique `celltype` entries.")
    }
    figure <- figure + scale_fill_manual(values = colors)
  }

  gginnards::delete_layers(x = figure, match_type = "GeomText")
}

