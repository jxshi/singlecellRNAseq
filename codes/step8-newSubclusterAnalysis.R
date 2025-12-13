# Step 8: new sub-cluster analysis for selected main cell types

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

base_dir <- "/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype"
setwd(base_dir)

assert_file <- function(path) {
  if (!file.exists(path)) stop(path, " not found")
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

save_marker_outputs <- function(obj, markers, prefix, assay = "RNA") {
  write.csv(markers, file = paste0(prefix, "_subcluster.markers.csv"))
  saveRDS(markers, file = paste0(prefix, "_subcluster.markers.Rds"))

  for (n in c(10, 6, 3)) {
    top_n_markers <- markers %>% group_by(cluster) %>% slice_max(n = n, order_by = avg_log2FC)
    heatmap <- DoHeatmap(obj, top_n_markers$gene, size = 3)
    ggsave(
      filename = paste0(prefix, "_DoHeatmap_top", n, ".pdf"),
      plot = heatmap, width = ifelse(n == 10, 14, 18), height = ifelse(n == 10, 20, 12)
    )

    dot <- DotPlot(obj, features = unique(top_n_markers$gene), assay = assay) + coord_flip()
    ggsave(
      plot = dot,
      filename = paste0(prefix, "_DotPlot_top", n, ".pdf"),
      device = cairo_pdf, width = 18, height = ifelse(n == 10, 24, 12)
    )
  }
}

save_dotplot <- function(object, genes, filename, assay = "RNA", width = 10, height = 8, title = NULL) {
  plot <- DotPlot(object, features = genes, assay = assay) + coord_flip()
  if (!is.null(title)) plot <- plot + ggtitle(title)
  ggsave(plot = plot, filename = filename, device = cairo_pdf, width = width, height = height)
  plot
}

run_harmony_workflow <- function(obj, npcs = 20, pre_resolution = 0.3, post_resolution = 0.5, k_param = 20) {
  obj <- obj %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA(npcs = npcs, verbose = FALSE)

  obj <- obj %>%
    RunUMAP(reduction = "pca", dims = 1:npcs, umap.method = "uwot") %>%
    RunTSNE(reduction = "pca", dims = 1:npcs) %>%
    FindNeighbors(reduction = "pca", k.param = k_param, dims = 1:npcs) %>%
    FindClusters(resolution = pre_resolution)

  pre_umap <- DimPlot(obj, reduction = "umap", split.by = "group", label = TRUE) +
    plot_annotation(title = "UMAP before harmony integration")
  pre_tsne <- DimPlot(obj, reduction = "tsne", split.by = "group", label = TRUE) +
    plot_annotation(title = "tSNE before harmony integration")
  pre_combined <- pre_umap + pre_tsne

  obj <- obj %>% RunHarmony(group.by.vars = "group", assay.use = "RNA", reduction = "pca")

  obj <- obj %>%
    RunUMAP(reduction = "harmony", dims = 1:npcs, verbose = FALSE) %>%
    RunTSNE(reduction = "harmony", dims = 1:npcs, verbose = FALSE) %>%
    FindNeighbors(reduction = "harmony", k.param = k_param, dims = 1:npcs) %>%
    FindClusters(resolution = post_resolution)

  post_umap <- DimPlot(obj, reduction = "umap", split.by = "group", label = TRUE) +
    plot_annotation(title = "UMAP after harmony integration")
  post_tsne <- DimPlot(obj, reduction = "tsne", split.by = "group", label = TRUE) +
    plot_annotation(title = "tSNE after harmony integration")
  post_combined <- post_umap + post_tsne

  list(object = obj, pre_plot = pre_combined, post_plot = post_combined)
}

process_targets <- function(indices, folder_names, base_dir) {
  for (idx in indices) {
    folder <- folder_names[idx]
    message("Processing ", folder)
    setwd(file.path(base_dir, folder))

    sce_path <- paste0(folder, ".Rds")
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

    saveRDS(sce_sub, file = paste0(folder, ".sce.sub.harmony.Rds"))

    markers <- FindAllMarkers(sce_sub, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
    save_marker_outputs(sce_sub, markers, prefix = paste0(folder, "_rna"))

    cluster_by_sample <- table(sce_sub@meta.data$seurat_clusters, sce_sub@meta.data$orig.ident)
    write.table(cluster_by_sample, paste0(folder, "_celltype_per_group.tsv"), quote = FALSE, sep = "\t", row.names = TRUE)

    prop_plot <- SubPropPlot(sce_sub, "batch") + coord_flip()
    ggsave(
      prop_plot,
      filename = paste0(folder, "_subgroup_celltype_proportion_plot.pdf"),
      device = cairo_pdf, width = 10, height = 3
    )

    setwd(base_dir)
  }
}

assert_file("celltype_file_folder.Rds")
folder_names <- readRDS("celltype_file_folder.Rds")

# Indices correspond to selected main cell types for further sub-clustering
process_targets(c(2, 3, 7, 8), folder_names, base_dir)

# Mono/Macro marker review panels (optional if MonoMacro folder exists)
monomacro_idx <- which(folder_names == "MonoMacro")
if (length(monomacro_idx) > 0) {
  mono_folder <- folder_names[monomacro_idx]
  mono_obj <- readRDS(file.path(base_dir, mono_folder, paste0(mono_folder, ".sce.sub.harmony.Rds")))
  DefaultAssay(mono_obj) <- "RNA"

  genes_mono_full <- c(
    "CD14", "PTPRC", "CD68",
    "IRF1", "IRF2BP2", "IFI27L2", "IFI30", "IFITM3", "IRF5", "IFI44L", "MX1", "IFIT1", "IFI44", "ISG15",
    "IFIT3", "IFIT2", "IFI6", "IFI16", "IFI35", "IFIH1", "IRF7", "IRF8", "IRF3", "IFITM2", "ISG20L2",
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
    "HLA-DRB1", "HLA-DRA", "HLA-DQB1", "HLA-DPB1", "HLA-DPA1", "CD74", "HLA-DMA", "HLA-DMB", "HLA-DRB5",
    "FCGR3A", "FCGR1A", "FCGR3B", "FCGR2A", "FCGR2B", "FCGR2C",
    "CD163", "CD36", "MARCO", "MRC1",
    "SIRPA", "SIGLEC10", "LILRB1", "LILRB2", "PDCD1", "SLAMF7",
    "MPP7", "MMP9", "MMP8",
    "S100A4", "S100A6", "S100A8", "S100A9", "S100A10", "S100A11", "S100A12",
    "GAS5", "GAS7", "VEGFA", "VEGFB", "IL1B", "IL18",
    "MERTK", "AXL", "TIMD4", "SLC40A1", "HMOX1",
    "ANXA1", "HMGB2", "NLRP3", "AIF1", "SIGLEC1", "RNASE2",
    "TPM2", "LILRA4", "LAMP3", "IDO1", "IDO2", "CD1E", "CD1C"
  )
  save_dotplot(
    mono_obj, genes_mono_full,
    filename = file.path(base_dir, mono_folder, paste0(mono_folder, "_DotPlot_interested_markers.pdf")),
    width = 16, height = 24
  )

  universal_mhc <- c(
    "CD14", "PTPRC", "CD68", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
    "HLA-DRB1", "HLA-DRA", "HLA-DQB1", "HLA-DPB1", "HLA-DPA1", "CD74", "HLA-DMA", "HLA-DMB", "HLA-DRB5"
  )
  save_dotplot(
    mono_obj, universal_mhc,
    filename = file.path(base_dir, mono_folder, paste0(mono_folder, "_DotPlot_universal_MHC.pdf")),
    width = 9, height = 10, title = "Universal and MHC markers"
  )

  interferon_genes <- c(
    "IRF1", "IRF2BP2", "IFI27L2", "IFI30", "IFITM3", "IRF5", "IFI44L", "MX1", "IFIT1", "IFI44",
    "ISG15", "IFIT3", "IFIT2", "IFI6", "IFI16", "IFI35", "IFIH1", "IRF7", "IRF8", "IRF3", "IFITM2", "ISG20L2"
  )
  save_dotplot(
    mono_obj, interferon_genes,
    filename = file.path(base_dir, mono_folder, paste0(mono_folder, "_DotPlot_interferon.pdf")),
    width = 9, height = 10, title = "Interferon genes"
  )

  receptors_mmp <- c(
    "FCGR3A", "FCGR1A", "FCGR3B", "FCGR2A", "FCGR2B", "FCGR2C",
    "CD163", "CD36", "MARCO", "MRC1",
    "SIRPA", "SIGLEC10", "LILRB1", "LILRB2", "PDCD1", "SLAMF7",
    "MPP7", "MMP9", "MMP8"
  )
  save_dotplot(
    mono_obj, receptors_mmp,
    filename = file.path(base_dir, mono_folder, paste0(mono_folder, "_DotPlot_receptors_MMP.pdf")),
    width = 9, height = 10,
    title = "Fc receptors, scavenger receptors, don't-eat-me and MMP genes"
  )

  s100_cytokine <- c(
    "S100A4", "S100A6", "S100A8", "S100A9", "S100A10", "S100A11", "S100A12",
    "GAS5", "GAS7", "VEGFA", "VEGFB", "IL1B", "IL18",
    "MERTK", "AXL", "TIMD4", "SLC40A1", "HMOX1",
    "ANXA1", "HMGB2", "NLRP3", "AIF1", "SIGLEC1", "RNASE2"
  )
  save_dotplot(
    mono_obj, s100_cytokine,
    filename = file.path(base_dir, mono_folder, paste0(mono_folder, "_DotPlot_S100_cytokine.pdf")),
    width = 9, height = 10,
    title = "S100 family, cytokines and other markers"
  )
}

# B/plasma cell marker overview
if (dir.exists(file.path(base_dir, "Plasmacells"))) {
  setwd(file.path(base_dir, "Plasmacells"))
  plasma_file <- "Plasmacellssubcluster_best_resolution.Rds"
  assert_file(plasma_file)
  plasma_obj <- readRDS(plasma_file)
  plasma_genes <- c(
    "MS4A1", "CD79B", "CD79A", "CD22", "IGHM", "IGHD", "IGLC3", "CD19", "CD38",
    "SDC1", "TNFRSF17", "SLAMF7", "GPRC5D", "CD24", "GAS6", "GAS7", "IGF1", "CD27"
  )
  save_dotplot(
    plasma_obj, plasma_genes, assay = "SCT",
    filename = "DotPlot_check_interested_markers_for_B_Plasma_cells.pdf",
    width = 10, height = 8
  )
  setwd(base_dir)
}
