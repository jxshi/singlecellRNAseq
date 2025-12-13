PropPlot <- function(object, groupBy, colors = NULL) {
  library(dplyr)
  library(ggplot2)

  # Step 1. get data
  plot_data <- object@meta.data %>%
    dplyr::select({{ groupBy }}, celltype) %>%
    dplyr::rename(group = {{ groupBy }})

  # Step 2. base plot
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
    legend.title  = "Seurat Cluster",
    ggtheme       = ggpubr::theme_pubclean()
  ) +
    theme(
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(color = "black", lineend = "round"),
      legend.position = "right",
      axis.text.x = element_text(size = 15, color = "black", family = "Arial", angle=45, hjust=1,vjust=1),
      axis.text.y = element_text(size = 15, color = "black", family = "Arial"),
      legend.text  = element_text(family = "Arial", size = 10, color = "black"),
      legend.title = element_text(family = "Arial", size = 13, color = "black")
    )

  # Step 2b. override colors manually if provided
  if (!is.null(colors)) {
    figure <- figure + scale_fill_manual(values = colors)
  }

  # Step 3. remove some elements (text labels from ggstatsplot)
  figure <- gginnards::delete_layers(x = figure, match_type = "GeomText")

  return(figure)
}

