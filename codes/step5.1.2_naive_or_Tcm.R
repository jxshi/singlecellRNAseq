#################### Step 5.1.2. naive_or_Tcm cells ####################

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
  "Seurat", "tidyverse", "reticulate", "ggplot2", "patchwork", "SingleR", "dplyr",
  "celldex", "RColorBrewer", "SingleCellExperiment"
))

source("~/software/functions/SubClusterPropPlot.R")
source("~/software/functions/DoHeatmapPlot.R")
source("~/software/functions/convertHumanGeneList.R")

assert_file <- function(path) {
  if (!file.exists(path)) stop(path, " not found")
}

save_dotplot <- function(object, genes, filename, assay = "RNA", width = 10, height = 10, title = NULL) {
  plot <- DotPlot(object, features = genes, assay = assay, cols = c("lightgrey", "blue", "red", "green", "cyan")) +
    coord_flip() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1.2, vjust = 1.1, size = 8),
      strip.text  = element_text(size = 8)
    )
  if (!is.null(title)) plot <- plot + ggtitle(title)
  ggsave(filename, plot = plot, device = cairo_pdf, width = width, height = height)
  plot
}

run_harmony_workflow <- function(obj_list, prefix, dims = 1:30, res_pre = 0.8, res_post = 0.5) {
  prepped <- lapply(obj_list, function(obj) {
    DefaultAssay(obj) <- "RNA"
    NormalizeData(obj) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 3000)
  })

  merged <- merge(prepped[[1]], y = prepped[-1])
  DefaultAssay(merged) <- "RNA"

  merged <- merged %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA(npcs = max(dims), verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = dims, umap.method = "uwot") %>%
    RunTSNE(reduction = "pca", dims = dims) %>%
    FindNeighbors(reduction = "pca", k.param = 20, dims = dims) %>%
    FindClusters(resolution = res_pre)

  before <- DimPlot(merged, reduction = "umap", split.by = "group", label = TRUE) +
    plot_annotation(title = "UMAP plot of samples before harmony integration")
  before_tsne <- DimPlot(merged, reduction = "tsne", split.by = "group", label = TRUE) +
    plot_annotation(title = "TSNE plot of samples before harmony integration")
  ggsave("umap_tsne_plots_before_using_harmony.pdf", plot = before + before_tsne, device = "pdf", width = 24, height = 12)

  merged <- merged %>%
    RunHarmony(group.by.vars = c("orig.ident", "group"), assay.use = "RNA", reduction = "pca") %>%
    RunUMAP(reduction = "harmony", dims = dims, verbose = FALSE) %>%
    RunTSNE(reduction = "harmony", dims = dims, verbose = FALSE) %>%
    FindNeighbors(reduction = "harmony", k.param = 20, dims = dims) %>%
    FindClusters(resolution = res_post)

  after <- DimPlot(merged, reduction = "umap", split.by = "group", label = TRUE) +
    plot_annotation(title = "UMAP plot of samples after harmony integration")
  after_tsne <- DimPlot(merged, reduction = "tsne", split.by = "group", label = TRUE) +
    plot_annotation(title = "TSNE plot of samples after harmony integration")
  ggsave("umap_tsne_plots_after_using_harmony.pdf", plot = after + after_tsne, device = "pdf", width = 24, height = 12)

  saveRDS(merged, file = paste0(prefix, ".sce.sub.harmony.Rds"))
  merged
}

export_markers_and_plots <- function(obj, prefix, assay = "RNA", min_pct = 0.25, logfc = 0.25, tops = c(10, 6, 3)) {
  DefaultAssay(obj) <- assay
  markers <- FindAllMarkers(obj, only.pos = FALSE, min.pct = min_pct, logfc.threshold = logfc)
  write.csv(markers, file = paste0(prefix, "_rna_subcluster.markers.csv"))
  saveRDS(markers, file = paste0(prefix, "_rna_subcluster.markers.Rds"))

  invisible(lapply(tops, function(n) {
    topn <- markers %>% group_by(cluster) %>% slice_max(n = n, order_by = avg_log2FC)
    ggsave(
      filename = paste0(prefix, "_rna_DoHeatmap_top", n, "_markers_by_clusters.pdf"),
      plot = DoHeatmap(obj, features = topn$gene, size = 3), width = 18, height = ifelse(n == 10, 20, 12)
    )
    ggsave(
      plot = DotPlot(obj, features = unique(topn$gene), assay = assay) + coord_flip(),
      filename = paste0(prefix, "_rna_DotPlot_check_top", n, "_markers_by_clusters.pdf"),
      device = cairo_pdf, width = 18, height = ifelse(n == 10, 24, 12)
    )
  }))

  markers
}

apply_cluster_labels <- function(obj, mapping, column) {
  obj[[column]] <- "NA"
  for (i in seq_len(nrow(mapping))) {
    obj@meta.data[obj$seurat_clusters == mapping$ClusterID[i], column] <- mapping[[column]][i]
  }
  obj
}

write_proportion_tables <- function(obj, column, prefix) {
  cltyPerIdent <- table(obj@meta.data[[column]], obj@meta.data$orig.ident)
  write.table(cltyPerIdent, paste0(prefix, "_", column, "_per_orig.ident.tsv"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

  cltyPerGroup <- table(obj@meta.data[[column]], obj@meta.data$group)
  write.table(cltyPerGroup, paste0(prefix, "_", column, "_per_group.tsv"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
}

working_dir <- "/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/3-celltype/Tcells/naiveTcm"
if (dir.exists(working_dir)) setwd(working_dir)

prefix <- "naiveTcm"
dataset_path <- "sce.naiveTcm.Rds"
assert_file(dataset_path)

sce <- readRDS(dataset_path)
DefaultAssay(sce) <- "RNA"

sce_list <- SplitObject(sce, split.by = "group")
sce.sub <- run_harmony_workflow(sce_list, prefix)

sce_markers <- export_markers_and_plots(sce.sub, prefix)

# CD4/CD8 panels
cd_markers <- convertHumanGeneList(c("CD4", "CD8A", "CD3D", "CD3E", "TRAC", "IL7R", "CD44", "IL2RA", "FOXP3", "CTLA4", "GZMB", "PRF1", "CCL5", "NKG7", "GZMA"))
save_dotplot(sce.sub, cd_markers, paste0(prefix, "_Dotplot_plot.pdf"), width = 10, height = 10)

gd_markers <- c("Tcrg-C1", "Tcrg-C2", "Tcrg-C3", "Tcrg-C4", "Trac", "Trbc1", "Trbc2")
save_dotplot(sce.sub, gd_markers, paste0(prefix, "_gdT_markers.pdf"), width = 5, height = 8)

# Annotate clusters
subcelltype <- data.frame(
  ClusterID = 0:10,
  subcelltype = c(
    "CD4_Tcm", "CD8_Tem", "CD8_Tem", "CD4_Tem", "Tem",
    "eTreg", "CD4_Tcm", "CD4_Tcm", "CD8_Tem", "CD8_Teff", "CD8_Tem"
  )
)

sce.sub <- apply_cluster_labels(sce.sub, subcelltype, "subcelltype")
write_proportion_tables(sce.sub, "subcelltype", prefix)
Idents(sce.sub) <- sce.sub$subcelltype
saveRDS(sce.sub, "naiveTcm_w_subcelltype.Rds")

# Plot proportion plot on subclusters
cltyPerSubroup <- table(sce.sub@meta.data$seurat_clusters, sce.sub@meta.data$orig.ident)
write.table(cltyPerSubroup, paste0(prefix, "_celltype_per_group.tsv"), quote = FALSE, sep = "\t", row.names = TRUE)

ggsave(
  filename = paste0(prefix, "_subgroup_celltype_proportion_plot.pdf"),
  plot = SubPropPlot(sce.sub, "batch") + coord_flip(),
  device = cairo_pdf, width = 10, height = 3
)

# Mono/Macro marker review panels
mono_macro_panels <- list(
  list(
    title = "Universal, MHC and innate receptors",
    genes = c(
      "CD14", "PTPRC", "CD68",
      "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
      "HLA-DRB1", "HLA-DRA", "HLA-DQB1", "HLA-DPB1", "HLA-DPA1", "CD74", "HLA-DMA", "HLA-DMB", "HLA-DRB5"
    ),
    filename = paste0(prefix, "_DotPlot_univ_mhc_by_clusters.pdf"), width = 9, height = 10
  ),
  list(
    title = "Interferon genes",
    genes = c(
      "IRF1", "IRF2BP2", "IFI27L2", "IFI30", "IFITM3", "IRF5", "IFI44L", "MX1", "IFIT1", "IFI44", "ISG15", "IFIT3", "IFIT2", "IFI6",
      "IFI16", "IFI35", "IFIH1", "IRF7", "IRF8", "IRF3", "IFITM2", "ISG20L2"
    ),
    filename = paste0(prefix, "_DotPlot_interferon_markers_by_clusters.pdf"), width = 9, height = 10
  ),
  list(
    title = "Fc/Scavenger receptors and enzymes",
    genes = c(
      "FCGR3A", "FCGR1A", "FCGR3B", "FCGR2A", "FCGR2B", "FCGR2C",
      "CD163", "CD36", "MARCO", "MRC1",
      "SIRPA", "SIGLEC10", "LILRB1", "LILRB2", "PDCD1", "SLAMF7",
      "MPP7", "MMP9", "MMP8"
    ),
    filename = paste0(prefix, "_DotPlot_fcReceptor_markers_by_clusters.pdf"), width = 9, height = 10
  ),
  list(
    title = "S100 family, cytokines and other markers",
    genes = c(
      "S100A4", "S100A6", "S100A8", "S100A9", "S100A10", "S100A11", "S100A12",
      "GAS5", "GAS7", "VEGFA", "VEGFB", "IL1B", "IL18",
      "MERTK", "AXL", "TIMD4", "SLC40A1", "HMOX1",
      "ANXA1", "HMGB2", "NLRP3", "AIF1", "SIGLEC1", "RNASE2"
    ),
    filename = "DotPlot_s100_cytokine_markers_by_clusters.pdf", width = 9, height = 10
  )
)

invisible(lapply(mono_macro_panels, function(panel) {
  save_dotplot(
    sce.sub, panel$genes, filename = panel$filename, assay = "RNA",
    width = panel$width, height = panel$height, title = panel$title
  )
}))

# Optional plasma/B cell check
if (file.exists("Plasmacells/Plasmacellssubcluster_best_resolution.Rds")) {
  plasma_sce <- readRDS("Plasmacells/Plasmacellssubcluster_best_resolution.Rds")
  plasma_genes <- c("MS4A1", "CD79B", "CD79A", "CD22", "IGHM", "IGHD", "IGLC3", "CD19", "CD38", "SDC1", "TNFRSF17", "SLAMF7", "GPRC5D", "CD24",
                    "GAS6", "GAS7", "IGF1", "CD27")
  save_dotplot(plasma_sce, plasma_genes, "Plasmacells/DotPlot_check_interested_markers_for_B_Plasma_cells.pdf", assay = "SCT", width = 10, height = 8)
}
