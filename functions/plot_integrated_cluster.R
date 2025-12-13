plot_integrated_clusters <- function(srat,
                                     cluster_col = "subcelltype",  # was: "seurat_clusters"
                                     group_col   = "group") {

  ## take an integrated Seurat object, plot distributions over group_col
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)

  meta <- srat@meta.data

  if (!cluster_col %in% colnames(meta)) {
    stop(sprintf("Column '%s' not found in meta.data.", cluster_col))
  }
  if (!group_col %in% colnames(meta)) {
    stop(sprintf("Column '%s' not found in meta.data.", group_col))
  }

  # make sure cluster is factor
  meta[[cluster_col]] <- as.factor(meta[[cluster_col]])

  # table: rows = subcelltype (or whatever cluster_col is), columns = group
  count_table <- table(meta[[cluster_col]], meta[[group_col]])
  count_mtx   <- as.data.frame.matrix(count_table)
  count_mtx$cluster <- rownames(count_mtx)

  melt_mtx    <- melt(count_mtx)
  melt_mtx$cluster <- as.factor(melt_mtx$cluster)

  # total cells per subcelltype
  cluster_size <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)

  # order clusters by size (largest on top in the plot)
  ordered_levels <- cluster_size$cluster[order(cluster_size$value, decreasing = TRUE)]
  cluster_size$cluster <- factor(cluster_size$cluster, levels = ordered_levels)
  melt_mtx$cluster     <- factor(melt_mtx$cluster,     levels = ordered_levels)

  colnames(melt_mtx)[2] <- "dataset"  # rename group_col to "dataset" for plotting

  p1 <- ggplot(cluster_size, aes(y = cluster, x = value)) +
    geom_bar(position = "dodge", stat = "identity", fill = "grey60") +
    theme_bw() +
    scale_x_log10() +
    xlab("Cells per subcelltype (log10 scale)") +
    ylab("")

  p2 <- ggplot(melt_mtx, aes(x = cluster, y = value, fill = dataset)) +
    geom_bar(position = "fill", stat = "identity") +
    theme_bw() +
    coord_flip() +
    ggsci::scale_fill_nejm() +
    ylab("Fraction of cells in each dataset") +
    xlab(cluster_col) +
    theme(legend.position = "top")

  p2 + p1 + plot_layout(widths = c(3, 1))
}

