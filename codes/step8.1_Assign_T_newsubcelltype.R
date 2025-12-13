# Step 8.1: assign refined T-cell subtypes

load_required_packages <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
  }
  invisible(lapply(pkgs, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
}

load_required_packages(c(
  "Seurat", "tidyverse", "ggplot2", "patchwork", "SingleR", "dplyr", "celldex",
  "RColorBrewer", "SingleCellExperiment", "harmony"
))

source("~/software/functions/SubClusterPropPlot.R")
source("~/software/functions/DoHeatmapPlot.R")
source("~/software/functions/convertHumanGeneList.R")
source("~/software/functions/custom_seurat_functions.R")
source("~/software/functions/groupbarcharts.R")
source("~/software/functions/groupbarchartssubcluster.R")
source("~/software/functions/PropPlot.R")

assert_file <- function(path) {
  if (!file.exists(path)) stop(path, " not found")
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
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
  write.csv(markers, file = paste0(prefix, "_subcluster.markers.csv"))
  saveRDS(markers, file = paste0(prefix, "_subcluster.markers.Rds"))

  for (n in c(10, 6, 3)) {
    top_markers <- markers %>% group_by(cluster) %>% slice_max(n = n, order_by = avg_log2FC)
    heatmap <- DoHeatmap(obj, top_markers$gene, size = 3)
    ggsave(paste0(prefix, "_DoHeatmap_top", n, ".pdf"), plot = heatmap, width = 18, height = ifelse(n == 10, 24, 12))

    dot <- DotPlot(obj, features = unique(top_markers$gene), assay = assay) + coord_flip()
    ggsave(paste0(prefix, "_DotPlot_top", n, ".pdf"), plot = dot, device = cairo_pdf, width = 18, height = ifelse(n == 10, 24, 12))
  }
}

save_dotplot <- function(object, genes, filename, assay = "RNA", width = 10, height = 8, title = NULL) {
  plot <- DotPlot(object, features = genes, assay = assay) + coord_flip()
  if (!is.null(title)) plot <- plot + ggtitle(title)
  ggsave(plot = plot, filename = filename, device = cairo_pdf, width = width, height = height)
  plot
}

run_harmony_workflow <- function(obj, npcs = 30, resolution = 1.5, k_param = 20) {
  obj <- obj %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA(npcs = npcs, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:npcs, umap.method = "uwot") %>%
    RunTSNE(reduction = "pca", dims = 1:npcs) %>%
    FindNeighbors(reduction = "pca", k.param = k_param, dims = 1:npcs) %>%
    FindClusters(resolution = resolution)

  pre_plots <- DimPlot(obj, reduction = "umap", split.by = "group", label = TRUE) +
    DimPlot(obj, reduction = "tsne", split.by = "group", label = TRUE)

  obj <- obj %>% RunHarmony(group.by.vars = "group", assay.use = "RNA", reduction = "pca")

  obj <- obj %>%
    RunUMAP(reduction = "harmony", dims = 1:npcs, verbose = FALSE) %>%
    RunTSNE(reduction = "harmony", dims = 1:npcs, verbose = FALSE) %>%
    FindNeighbors(reduction = "harmony", k.param = k_param, dims = 1:npcs) %>%
    FindClusters(resolution = resolution)

  post_plots <- DimPlot(obj, reduction = "umap", split.by = "group", label = TRUE) +
    DimPlot(obj, reduction = "tsne", split.by = "group", label = TRUE)

  list(object = obj, pre_plot = pre_plots, post_plot = post_plots)
}

base_dir <- "/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype"
setwd(base_dir)
assert_file("celltype_file_folder.Rds")
folder_names <- readRDS("celltype_file_folder.Rds")

t_folder <- folder_names[1]
setwd(file.path(base_dir, t_folder))

sce_path <- paste0(t_folder, ".Rds")
assert_file(sce_path)
sce <- readRDS(sce_path)
DefaultAssay(sce) <- "RNA"

sce_list <- SplitObject(sce, split.by = "group")
sce_list <- lapply(sce_list, function(obj) {
  DefaultAssay(obj) <- "RNA"
  NormalizeData(obj) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 3000)
})

sce_sub <- Reduce(function(x, y) merge(x, y), sce_list)
DefaultAssay(sce_sub) <- "RNA"

harmony_res <- run_harmony_workflow(sce_sub)
sce_sub <- harmony_res$object
ggsave("umap_tsne_plots_before_using_harmony.pdf", plot = harmony_res$pre_plot, device = "pdf", width = 24, height = 12)
ggsave("umap_tsne_plots_after_using_harmony.pdf", plot = harmony_res$post_plot, device = "pdf", width = 24, height = 12)

saveRDS(sce_sub, file = paste0(t_folder, ".sce.sub.harmony.Rds"))

markers <- FindAllMarkers(sce_sub, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
save_marker_outputs(sce_sub, markers, prefix = paste0(t_folder, "_rna"))

cluster_by_sample <- table(sce_sub@meta.data$seurat_clusters, sce_sub@meta.data$orig.ident)
write.table(cluster_by_sample, paste0(t_folder, "_celltype_per_group.tsv"), quote = FALSE, sep = "\t", row.names = TRUE)

prop_plot <- SubPropPlot(sce_sub, "batch") + coord_flip()
ggsave(
  prop_plot,
  filename = paste0(t_folder, "_subgroup_celltype_proportion_plot.pdf"),
  device = cairo_pdf, width = 10, height = 3
)

newsub_labels <- data.frame(
  ClusterID = 0:26,
  newsubcelltype = c(
    "Tn", "Cd8 Teff", "Cd8 Tpex", "Th1", "Treg", "Cd8 Teff", "Cd8 Tm", "Tfh", "NK",
    "Tn", "Cd4 T early activated", "Proliferating Cd8 Teff", "Proliferating Cd8 T",
    "Cd4 Tem", "gdT", "Cd8 Tex", "Cd8 Teff", "Th2", "Proliferating Cd8 T", "Treg",
    "IFN-responsive Cd8 Teff", "Activated Th2", "Proliferating Treg", "Ly6c+ Cd8 Teff",
    "GC B", "GC B", "Cd4 Tcm"
  )
)

main_labels <- data.frame(
  ClusterID = 0:26,
  maincelltype = c(
    rep("T cells", 8), "NK cells", rep("T cells", 15), rep("Unknown", 2), "T cells"
  )
)

sce_sub <- apply_cluster_labels(sce_sub, newsub_labels, "newsubcelltype")
write_proportion_tables(sce_sub, "newsubcelltype", prefix = "newsubcelltype")
sce_sub <- apply_cluster_labels(sce_sub, main_labels, "maincelltype")
write_proportion_tables(sce_sub, "maincelltype", prefix = "maincelltype")

sce_t <- sce_sub
saveRDS(sce_t, file.path(base_dir, "sce.T_w_correction.Rds"))

combined_colors <- c(
  "#3C8A32", "#28357B", "#0F7374", "#116D69", "#5EACA7", "#B8D8DF", "#D7E9EC",
  "#9BC9D9", "#BFC7E5", "#A9B4D7", "#D8A5CE", "#BE84B8", "#87387D", "#D2A21B",
  "#F5C34D", "#C0261D", "#A62071", "#F1A800", "#DDA873", "#1A1B4B", "#FF5733"
)

# Plot figures
graphs_dir <- file.path(base_dir, "graphs", "TcellsFigures")
ensure_dir(graphs_dir)
setwd(graphs_dir)

Idents(sce_t) <- sce_t$newsubcelltype
sce_plot <- subset(sce_t, idents = setdiff(levels(Idents(sce_t)), c("GC B", "NK")))
Idents(sce_plot) <- sce_plot$newsubcelltype
sce_plot$subcelltype <- sce_plot$newsubcelltype
sce_plot$group <- factor(sce_plot$group, levels = c("PBS", "DD_mGE", "DD_mIL12", "DD_mGE12", "TD_mGE12"))

sce_plot$subcelltype <- factor(
  sce_plot$subcelltype,
  levels = c(
    "Tn", "Cd4 Tcm", "Cd4 Tem", "Th1", "Th2", "Tfh", "Treg", "Proliferating Treg",
    "Cd8 Tm", "Cd8 Tem", "Proliferating Cd8 T", "Proliferating Cd8 Teff", "Cd8 Teff",
    "IFN-responsive Cd8 Teff", "Ly6c+ Cd8 Teff", "Cd8 Tex", "gdT"
  )
)

p <- DimPlot(sce_plot, group.by = "subcelltype", raster = FALSE, label = TRUE, cols = combined_colors)
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
markers_plot <- FindAllMarkers(sce_plot, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
save_marker_outputs(sce_plot, markers_plot, prefix = "cca")

# Marker review panels
save_dotplot(
  sce_plot,
  c("Cd3d", "Cd3e", "Cd4", "Cd8a", "Cd44", "Il7r", "Sell", "Lef1", "Ccr7", "Cxcr5", "Tcf7", "Gzmk", "Il17a", "Gata3", "Il4"),
  "DotPlot_core_T_markers.pdf"
)

save_dotplot(
  sce_plot,
  c("Trgc1", "Trgc2", "Trdc", "Trgv1", "Trgv2", "Ncr1", "Klrb1c", "Itgam", "Klrc1"),
  "DotPlot_gdT_NK_markers.pdf"
)

exhausted_markers <- c(
  "Pdcd1", "Lag3", "Havcr2", "Ctla4", "Tigit", "Cd244", "Cd160", "Entpd1",
  "Tox", "Tox2", "Nr4a1", "Nr4a2", "Eomes", "Tcf7",
  "Cxcr6", "Cxcl13", "Tnfrsf9", "Tnfrsf18", "Icos", "Cd27", "Cd28", "Cxcr5", "Il7r"
)
save_dotplot(sce_plot, exhausted_markers, "DotPlot_exhaustion_markers.pdf")

pre_exhausted <- c(
  "Tcf7", "Slamf6", "Cxcr5", "Il7r", "Pdcd1", "Tox", "Cd28", "Bcl6", "Btg1", "Btg2",
  "Ccr7", "Sell", "Lef1", "Id3", "Eomes", "Tnfrsf4", "Tnfrsf9", "Icos"
)
save_dotplot(sce_plot, pre_exhausted, "DotPlot_pre_exhausted.pdf")

cm_markers <- convertHumanGeneList(c("IL7R", "CCR7", "SELL", "TCF7", "LEF1", "BCL2", "GZMK", "CXCR5", "IL17"))
save_dotplot(sce_plot, cm_markers, "DotPlot_central_memory_markers.pdf")

effector_cd8 <- c(
  "Gzmb", "Gzma", "Gzmk", "Prf1", "Ifng", "Tnf", "Ccl3", "Ccl4", "Ccl5", "Nkg7", "Ctsw", "Cst7",
  "Klrc1", "Klrc2", "Klrd1", "Klrk1", "Hopx", "Tbx21", "Zeb2"
)
save_dotplot(sce_plot, effector_cd8, "DotPlot_effector_cd8.pdf")

nk_markers <- c(
  "Ncr1", "Klrk1", "Klrd1", "Klrc1", "Klrc2", "Klre1", "Tyrobp", "Fcer1g", "Fcgr3",
  "Gzmb", "Gzma", "Gzmk", "Prf1", "Ctsw", "Cst7", "Ifng", "Fasl", "Tnfsf10", "Xcl1", "Ccl3", "Ccl4", "Ccl5",
  "Tbx21", "Eomes", "Bhlhe40", "Irf8", "Il2rb", "Cd160", "Klrg1", "Cd69", "Serpinb9", "Hilpda", "Ddit4", "Hif1a", "Mki67"
)
save_dotplot(sce_plot, nk_markers, "DotPlot_NK_activation.pdf")

nk_unique <- c(
  "Ncr1", "Klrb1c", "Klra3", "Klra7", "Klra8", "Klra9", "Klra12", "Klra17", "Klrd1", "Klrc1", "Klrc2", "Klre1",
  "Tyrobp", "Fcer1g", "Fcgr3", "Gzmc", "Gzmg", "Gzme", "Eomes", "Irf8", "Cd160", "Cd244"
)
save_dotplot(sce_plot, nk_unique, "DotPlot_NK_unique.pdf")

treg_markers <- c(
  "Foxp3", "Il2ra", "Ctla4", "Ikzf2", "Tnfrsf18", "Tnfrsf4", "Ccr7", "Ccr4", "Ccr8", "Cxcr3",
  "Tigit", "Havcr2", "Lag3", "Pdcd1", "Entpd1", "Nt5e", "Il10", "Tgfb1", "Ebi3", "Il12a", "Nrp1", "Bcl2", "Bcl2l1",
  "Gzmb", "Layn", "Il1rl1", "Areg"
)
save_dotplot(sce_plot, treg_markers, "DotPlot_Treg.pdf")

cd8_exhausted <- c(
  "Cd8a", "Cd8b1", "Pdcd1", "Lag3", "Havcr2", "Cd160", "Entpd1", "Tnfrsf9",
  "Gzmf", "Gzmd", "Gzmc", "Gzmb", "Prf1", "Nkg7", "Cst7", "Klrd1", "Klrc1", "Klrc2", "Klrk1",
  "Ifng", "Tbx21", "S100a4", "Bhlhe40", "Nr4a2"
)
save_dotplot(sce_plot, cd8_exhausted, "DotPlot_cd8_exhausted.pdf")

myeloid_mixed <- c(
  "S100a8", "S100a9", "Ngp", "Retnlg", "Cxcl1", "Cxcl2", "Lyz2", "Cd14", "C1qa", "Aif1", "Lst1", "Apoe", "Cebpb",
  "Alox5ap", "Fcer1g", "Tyrobp", "Ctss", "Ctsd", "Acp5", "H2-Ab1", "H2-Aa", "H2-Eb1", "Cd74",
  "Cd3e", "Cd8a", "Cd8b1", "Gzmk", "Nkg7", "Klrc1", "Klrc2", "Ccl5", "Cxcr3"
)
save_dotplot(sce_plot, myeloid_mixed, "DotPlot_myeloid_mixed.pdf")

th1_markers <- c("Tbx21", "Ifng", "Il2", "Cxcr3", "Ccr5", "Stat1", "Stat4", "Irf1", "Irf8", "Ccl3", "Ccl4", "Ccl5", "Gzmb", "Nkg7")
th2_markers <- c("Gata3", "Il4", "Il5", "Il13", "Maf", "Il1rl1", "Il7r", "Ccr4", "Ccl1", "Ccl17", "Ccl22", "Areg", "Csf2", "Ramp1", "Ramp3")
save_dotplot(sce_plot, unique(c(th1_markers, th2_markers)), "DotPlot_Th1_Th2.pdf")

setwd(base_dir)
