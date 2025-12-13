# Step 8.2: assign monocyte/macrophage subtypes

load_required_packages <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
  }
  invisible(lapply(pkgs, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
}

load_required_packages(c("Seurat", "dplyr", "ggplot2", "ggsci"))

source("~/software/functions/groupbarcharts.R")
source("~/software/functions/groupbarchartssubcluster.R")
source("~/software/functions/PropPlot.R")
source("~/software/functions/custom_seurat_functions.R")

assert_file <- function(path) {
  if (!file.exists(path)) stop(path, " not found")
}

apply_cluster_labels <- function(obj, mapping, column) {
  obj[[column]] <- "NA"
  for (i in seq_len(nrow(mapping))) {
    obj@meta.data[obj$seurat_clusters == mapping$ClusterID[i], column] <- mapping[[column]][i]
  }
  obj
}

write_proportion_tables <- function(obj, column, prefix = column) {
  per_sample <- table(obj@meta.data[[column]], obj@meta.data$orig.ident)
  write.table(per_sample, paste0(prefix, "_per_orig.ident.tsv"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

  per_group <- table(obj@meta.data[[column]], obj@meta.data$group)
  write.table(per_group, paste0(prefix, "_per_group.tsv"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
}

save_marker_outputs <- function(obj, markers, prefix, assay = "RNA") {
  write.csv(markers, file = paste0(prefix, "_markers.csv"))
  saveRDS(markers, file = paste0(prefix, "_markers.Rds"))

  for (n in c(10, 6, 3)) {
    top_markers <- markers %>% group_by(cluster) %>% slice_max(n = n, order_by = avg_log2FC)
    heatmap <- DoHeatmap(obj, top_markers$gene, size = 3)
    ggsave(paste0(prefix, "_DoHeatmap_top", n, ".pdf"), plot = heatmap, width = 18, height = ifelse(n == 10, 24, 12))

    dot <- DotPlot(obj, features = unique(top_markers$gene), assay = assay) + coord_flip()
    ggsave(paste0(prefix, "_DotPlot_top", n, ".pdf"), plot = dot, device = "pdf", width = 18, height = ifelse(n == 10, 24, 12))
  }
}

save_dotplot <- function(object, genes, filename, assay = "RNA", width = 10, height = 8, title = NULL) {
  plot <- DotPlot(object, features = genes, assay = assay) + coord_flip()
  if (!is.null(title)) plot <- plot + ggtitle(title)
  ggsave(plot = plot, filename = filename, device = cairo_pdf, width = width, height = height)
  plot
}

base_dir <- "/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype"
mono_dir <- file.path(base_dir, "MonoMacro")
setwd(mono_dir)

sce_path <- file.path(mono_dir, "MonoMacro.sce.sub.harmony.Rds")
assert_file(sce_path)
sce_sub <- readRDS(sce_path)
sce_sub$group <- factor(sce_sub$group, levels = c("PBS", "DD_mGE", "DD_mIL12", "DD_mGE12", "TD_mGE12"))

DefaultAssay(sce_sub) <- "RNA"
Idents(sce_sub) <- sce_sub$seurat_clusters

m0_genes <- c("Csf1r", "Cd68", "Aif1", "Lyz2", "Adgre1", "Mertk", "Tyrobp", "C1qa", "C1qb", "C1qc", "Lgals3", "Cd36", "Ccl5", "H2-Aa", "H2-Ab1")
m1_genes <- c("Nos2", "Il1b", "Tnf", "Il6", "Cxcl9", "Cxcl10", "Stat1", "Irf5", "Cd86", "Cd80", "Gbp2", "Gbp5", "Socs3", "Ifit1", "Ifit3", "Isg15", "Rsad2", "Irf7")
m2_genes <- c("Mrc1", "Cd163", "Arg1", "Chil3", "Retnla", "Clec10a", "Il10", "Tgfb1", "Mgl2", "Ccl17", "Ccl22", "Vegfa", "Mmp12", "Spp1", "Apoe", "Pdgfb", "Timp2", "Fn1", "Anxa1", "Stab1", "Vsig4", "Klf4", "Stat6", "Irf4", "Cd63", "Il10ra", "Il4ra", "Ly6c1", "Ly6c2", "Slfn4", "Cd274", "Acp5", "Ctsk")

save_dotplot(
  sce_sub, genes = list(M0 = m0_genes, M1 = m1_genes, M2 = m2_genes),
  filename = "DotPlot_macrophage_states.pdf", title = "M0/M1/M2 marker panels"
)

sce_sub <- AddModuleScore(sce_sub, features = list(M0 = m0_genes), name = "M0_Score")
sce_sub <- AddModuleScore(sce_sub, features = list(M1 = m1_genes), name = "M1_Score")
sce_sub <- AddModuleScore(sce_sub, features = list(M2 = m2_genes), name = "M2_Score")

FeaturePlot(sce_sub, features = c("M0_Score1", "M1_Score1", "M2_Score1"))
VlnPlot(sce_sub, features = c("M0_Score1", "M1_Score1", "M2_Score1"), group.by = "seurat_clusters")

new_labels <- data.frame(
  ClusterID = 0:9,
  newsubcelltype = c(
    "Macrophages0", "Macrophages1", "Macrophages2", "Macrophages3", "Macrophages4",
    "Macrophages5", "Macrophages6", "Macrophages7", "Macrophages8", "Macrophages9"
  )
)

main_labels <- data.frame(ClusterID = 0:9, maincelltype = rep("Macrophages", 10))

sce_sub <- apply_cluster_labels(sce_sub, new_labels, "newsubcelltype")
write_proportion_tables(sce_sub, "newsubcelltype", prefix = "newsubcelltype")
sce_sub <- apply_cluster_labels(sce_sub, main_labels, "maincelltype")
write_proportion_tables(sce_sub, "maincelltype", prefix = "maincelltype")

sce_monomacro <- sce_sub
saveRDS(sce_monomacro, file.path(base_dir, "sce.MonoMacro_w_correction.Rds"))

combined_colors <- c(
  "#BC3C29", "#0072B5", "#E18727", "#008B45", "#5F559B", "#6F99AD", "#A20056", "#EE4C97", "#3B4992", "#008280", "#631879", "#C70039", "#FFDC91"
)

fig_dir <- file.path(base_dir, "graphs", "MonoMacroFigures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
setwd(fig_dir)

sce_plot <- sce_monomacro
Idents(sce_plot) <- sce_plot$newsubcelltype
sce_plot$subcelltype <- sce_plot$newsubcelltype

p <- DimPlot(sce_plot, group.by = "maincelltype", raster = FALSE, label = TRUE, cols = combined_colors)
ggsave("main_celltype_UMAP.pdf", plot = p, device = cairo_pdf, width = 8, height = 8)

p <- DimPlot(sce_plot, group.by = "group", raster = FALSE, cols = combined_colors)
ggsave("UMAP_by_group.pdf", plot = p, device = cairo_pdf, width = 8, height = 8)

p <- groupbarcharts(sce_plot, "group")
ggsave("compare_maincelltype_percentage_of_each_group.pdf", plot = p, device = cairo_pdf, width = 8, height = 8)

tmp <- table(sce_plot$group, sce_plot$subcelltype)
write.table(tmp, file = "subcelltype_per_group.tsv", quote = FALSE, col.names = NA, sep = "\t")

p <- groupbarchartssubcluster(sce_plot, "group")
ggsave("compare_percentage_of_cells_in_each_group.pdf", plot = p, device = cairo_pdf, width = 14, height = 8)

p <- plot_integrated_clusters(sce_plot, cluster_col = "subcelltype", group_col = "group")
ggsave("compare_percentage_of_subcelltype.pdf", plot = p, device = cairo_pdf, width = 8, height = 8)

Idents(sce_plot) <- sce_plot$subcelltype
DefaultAssay(sce_plot) <- "RNA"
sce_markers <- FindAllMarkers(sce_plot, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
save_marker_outputs(sce_plot, sce_markers, prefix = "cca")

setwd(base_dir)
