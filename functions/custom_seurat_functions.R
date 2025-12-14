plot_integrated_clusters <- function (srat) {
  if (!all(c("seurat_clusters", "group") %in% colnames(srat@meta.data))) {
    stop("`srat@meta.data` must contain `seurat_clusters` and `group` columns.")
  }

  count_table <- table(srat@meta.data$seurat_clusters, srat@meta.data$group)
  count_mtx   <- as.data.frame.matrix(count_table)
  count_mtx$cluster <- rownames(count_mtx)
  melt_mtx    <- reshape2::melt(count_mtx)
  melt_mtx$cluster <- as.factor(melt_mtx$cluster)

  cluster_size   <- stats::aggregate(value ~ cluster, data = melt_mtx, FUN = sum)

  sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)), decreasing = TRUE))
  cluster_size$cluster <- factor(cluster_size$cluster, levels = sorted_labels)
  melt_mtx$cluster <- factor(melt_mtx$cluster, levels = sorted_labels)

  colnames(melt_mtx)[2] <- "dataset"

  p1 <- ggplot2::ggplot(cluster_size, ggplot2::aes(y = cluster, x = value)) +
    ggplot2::geom_bar(position = "dodge", stat = "identity", fill = "grey60") +
    ggplot2::theme_bw() + ggplot2::scale_x_log10() +
    ggplot2::xlab("Cells per cluster, log10 scale") + ggplot2::ylab("")

  p2 <- ggplot2::ggplot(melt_mtx, ggplot2::aes(x = cluster, y = value, fill = dataset)) +
    ggplot2::geom_bar(position = "fill", stat = "identity") +
    ggplot2::theme_bw() + ggplot2::coord_flip() +
    ggsci::scale_fill_nejm() +
    ggplot2::ylab("Fraction of cells in each dataset") +
    ggplot2::xlab("Cluster number") +
    ggplot2::theme(legend.position = "top")

  p2 + p1 + patchwork::plot_layout(widths = c(3, 1))
}

plot_integrated_subcelltype <- function (srat) {
  if (!all(c("subcelltype", "group") %in% colnames(srat@meta.data))) {
    stop("`srat@meta.data` must contain `subcelltype` and `group` columns.")
  }

  count_table <- table(srat@meta.data$subcelltype, srat@meta.data$group)
  count_mtx   <- as.data.frame.matrix(count_table)
  count_mtx$cluster <- rownames(count_mtx)
  melt_mtx    <- reshape2::melt(count_mtx)
  melt_mtx$cluster <- as.factor(melt_mtx$cluster)

  cluster_size   <- stats::aggregate(value ~ cluster, data = melt_mtx, FUN = sum)

  sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)), decreasing = TRUE))
  cluster_size$cluster <- factor(cluster_size$cluster, levels = sorted_labels)
  melt_mtx$cluster <- factor(melt_mtx$cluster, levels = sorted_labels)

  colnames(melt_mtx)[2] <- "dataset"

  p1 <- ggplot2::ggplot(cluster_size, ggplot2::aes(y = cluster, x = value)) +
    ggplot2::geom_bar(position = "dodge", stat = "identity", fill = "grey60") +
    ggplot2::theme_bw() + ggplot2::scale_x_log10() +
    ggplot2::xlab("Cells per celltype, log10 scale") + ggplot2::ylab("")

  p2 <- ggplot2::ggplot(melt_mtx, ggplot2::aes(x = cluster, y = value, fill = dataset)) +
    ggplot2::geom_bar(position = "fill", stat = "identity") +
    ggplot2::theme_bw() + ggplot2::coord_flip() +
    ggsci::scale_fill_nejm() +
    ggplot2::ylab("Fraction of cells in each dataset") + ggplot2::xlab("cell type") +
    ggplot2::theme(legend.position = "top")

  p2 + p1 + patchwork::plot_layout(widths = c(3, 1))
}

