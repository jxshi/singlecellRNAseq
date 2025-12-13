#################### Step 4. Sub-cluster Analysis ####################
silent_load <- function(pkgs) {
  invisible(lapply(pkgs, function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
      stop("Package not installed: ", pkg)
    }
  }))
}

save_plot <- function(plot_obj, filename, width = 12, height = 8, device = cairo_pdf) {
  ggsave(plot = plot_obj, filename = filename, device = device, width = width, height = height)
}

with_dir <- function(path, code) {
  old_dir <- setwd(path)
  on.exit(setwd(old_dir))
  force(code)
}

plot_dim_splits <- function(object, reduction = "umap", split_by = "group", title_prefix = "") {
  DimPlot(object, reduction = reduction, split.by = split_by, label = TRUE) +
    plot_annotation(title = sprintf("%s plot of samples %s integration", toupper(reduction), title_prefix))
}

save_dotplot <- function(object, features, filename, assay = "RNA", title = NULL,
                         width = 18, height = 12, flip = TRUE) {
  plot <- DotPlot(object, features = features, assay = assay)
  if (flip) {
    plot <- plot + coord_flip()
  }
  if (!is.null(title)) {
    plot <- plot + ggtitle(title)
  }
  save_plot(plot, filename, width = width, height = height)
}

save_heatmap <- function(object, genes, filename, width = 18, height = 12, size = 3) {
  plot <- DoHeatmap(object, features = genes, size = size)
  save_plot(plot, filename, width = width, height = height)
}

summarize_markers <- function(markers, top_n = 10) {
  markers %>% group_by(cluster) %>% top_n(top_n, avg_log2FC)
}

run_marker_checks <- function(sce_sub, markers, prefix, assay = "RNA", pro = "_rna") {
  marker_file_prefix <- paste0(prefix, pro)

  write.csv(markers, file = paste0(marker_file_prefix, "_subcluster.markers.csv"))
  saveRDS(markers, file = paste0(marker_file_prefix, "_subcluster.markers.Rds"))

  top_counts <- c(`10` = 24, `6` = 12, `3` = 12)
  for (n in names(top_counts)) {
    top_genes <- summarize_markers(markers, as.numeric(n))
    genes <- unique(top_genes$gene)

    save_heatmap(
      sce_sub,
      genes,
      filename = paste0(marker_file_prefix, "_DoHeatmap_check_top", n, "_markers_by_clusters.pdf"),
      width = 18,
      height = ifelse(as.numeric(n) == 10, 20, 12)
    )

    save_dotplot(
      sce_sub,
      features = genes,
      assay = assay,
      filename = paste0(marker_file_prefix, "DotPlot_check_top", n, "_markers_by_clusters.pdf"),
      width = 18,
      height = top_counts[[n]]
    )
  }
}

silent_load(c(
  "Seurat", "tidyverse", "ggplot2", "patchwork", "SingleR",
  "dplyr", "celldex", "RColorBrewer", "SingleCellExperiment", "reticulate"
))

source("~/software/functions/SubClusterPropPlot.R")
source("~/software/functions/DoHeatmapPlot.R")

if (!exists("cname") && file.exists("celltype_name.Rds")) {
  cname <- readRDS("celltype_name.Rds")
}
if (!exists("fname") && file.exists("celltype_file_folder.Rds")) {
  fname <- readRDS("celltype_file_folder.Rds")
}

if (!exists("fname")) {
  stop("`fname` must be available or celltype_file_folder.Rds must exist in the working directory.")
}

target_indices <- if (exists("target_indices")) target_indices else c(8)

process_cluster <- function(idx) {
  with_dir(fname[idx], {
    sce <- readRDS(paste0(fname[idx], ".Rds"))
    DefaultAssay(sce) <- "RNA"

    sce_list <- SplitObject(sce, split.by = "group")
    sce_list <- lapply(sce_list, function(obj) {
      DefaultAssay(obj) <- "RNA"
      obj %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst", nfeatures = 3000)
    })

    sce_sub <- merge(sce_list[[1]], y = sce_list[-1])
    DefaultAssay(sce_sub) <- "RNA"

    sce_sub <- sce_sub %>%
      NormalizeData() %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
      ScaleData() %>%
      RunPCA(npcs = 10, verbose = FALSE) %>%
      RunUMAP(reduction = "pca", dims = 1:10, umap.method = "uwot") %>%
      RunTSNE(reduction = "pca", dims = 1:10) %>%
      FindNeighbors(reduction = "pca", k.param = 20, dims = 1:10) %>%
      FindClusters(resolution = 0.3)

    before_umap <- plot_dim_splits(sce_sub, reduction = "umap", title_prefix = "before harmony")
    before_tsne <- plot_dim_splits(sce_sub, reduction = "tsne", title_prefix = "before harmony")
    save_plot(before_umap + before_tsne, "umap_tsne_plots_before_using_harmony.pdf", width = 24, height = 12)

    sce_sub <- sce_sub %>%
      RunHarmony(group.by.vars = c("orig.ident", "group"), assay.use = "RNA", reduction = "pca") %>%
      RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
      RunTSNE(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
      FindNeighbors(reduction = "harmony", k.param = 20, dims = 1:30) %>%
      FindClusters(resolution = 0.5)

    after_umap <- plot_dim_splits(sce_sub, reduction = "umap", title_prefix = "after harmony")
    after_tsne <- plot_dim_splits(sce_sub, reduction = "tsne", title_prefix = "after harmony")
    save_plot(after_umap + after_tsne, "umap_tsne_plots_after_using_harmony.pdf", width = 24, height = 12)

    saveRDS(sce_sub, file = paste0(fname[idx], ".sce.sub.harmony.Rds"))

    markers <- FindAllMarkers(sce_sub, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
    run_marker_checks(sce_sub, markers, prefix = fname[idx], assay = "RNA")
  })
}

invisible(lapply(target_indices, process_cluster))

# Plot propotion plot on subclusters.
for (i in seq_along(fname)) {
  with_dir(fname[i], {
    sce <- readRDS(paste0(fname[i], ".sce.sub.harmony.Rds"))
    cltyPerSubroup <- table(sce@meta.data$seurat_cluster, sce@meta.data$orig.ident)
    write.table(cltyPerSubroup, paste0(fname[i], "celltype_per_group.tsv"), quote = FALSE, sep = "\t", row.names = TRUE)

    prop_plot <- SubPropPlot(sce, "batch") + coord_flip()
    save_plot(prop_plot, filename = paste0(fname[i], "_subgroup_celltype_proportion_plot.pdf"), width = 10, height = 3)
  })
}

# For MonoMacro subcluster
plot_monocyte_markers <- function(sce_obj, prefix, assay = "RNA", pro = "_rna") {
  groups <- list(
    interested = c(
      'CD14', 'PTPRC', 'CD68',
      'IRF1', 'IRF2BP2', 'IFI27L2', 'IFI30', 'IFITM3', 'IRF5', 'IFI44L', 'MX1', 'IFIT1', 'IFI44', 'ISG15', 'IFIT3', 'IFIT2', 'IFI6',
      'IFI16', 'IFI35', 'IFIH1', 'IRF7', 'IRF8', 'IRF3', 'IFITM2', 'ISG20L2',
      'HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F',
      'HLA-DRB1', 'HLA-DRA', 'HLA-DQB1', 'HLA-DPB1', 'HLA-DPA1', 'CD74', 'HLA-DMA', 'HLA-DMB', 'HLA-DRB5',
      'FCGR3A', 'FCGR1A', 'FCGR3B', 'FCGR2A', 'FCGR2B', 'FCGR2C',
      'CD163', 'CD36', 'MARCO', 'MRC1',
      'SIRPA', 'SIGLEC10', 'LILRB1', 'LILRB2', 'PDCD1', 'SLAMF7',
      'MPP7', 'MMP9', 'MMP8',
      'S100A4', 'S100A6', 'S100A8', 'S100A9', 'S100A10', 'S100A11', 'S100A12',
      'GAS5', 'GAS7', 'VEGFA', 'VEGFB', 'IL1B', 'IL18',
      'MERTK', 'AXL', 'TIMD4', 'SLC40A1', 'HMOX1',
      'ANXA1', 'HMGB2', 'NLRP3', 'AIF1', 'SIGLEC1', 'RNASE2',
      'TPM2', 'LILRA4',
      'LAMP3', 'IDO1', 'IDO2',
      'CD1E', 'CD1C'
    ),
    univ_mhc = c(
      'CD14', 'PTPRC', 'CD68',
      'HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F',
      'HLA-DRB1', 'HLA-DRA', 'HLA-DQB1', 'HLA-DPB1', 'HLA-DPA1', 'CD74', 'HLA-DMA', 'HLA-DMB', 'HLA-DRB5'
    ),
    interferon = c(
      'IRF1', 'IRF2BP2', 'IFI27L2', 'IFI30', 'IFITM3', 'IRF5', 'IFI44L', 'MX1', 'IFIT1', 'IFI44', 'ISG15', 'IFIT3', 'IFIT2', 'IFI6',
      'IFI16', 'IFI35', 'IFIH1', 'IRF7', 'IRF8', 'IRF3', 'IFITM2', 'ISG20L2'
    ),
    receptors = c(
      'FCGR3A', 'FCGR1A', 'FCGR3B', 'FCGR2A', 'FCGR2B', 'FCGR2C',
      'CD163', 'CD36', 'MARCO', 'MRC1',
      'SIRPA', 'SIGLEC10', 'LILRB1', 'LILRB2', 'PDCD1', 'SLAMF7',
      'MPP7', 'MMP9', 'MMP8'
    ),
    s100_cytokine = c(
      'S100A4', 'S100A6', 'S100A8', 'S100A9', 'S100A10', 'S100A11', 'S100A12',
      'GAS5', 'GAS7', 'VEGFA', 'VEGFB', 'IL1B', 'IL18',
      'MERTK', 'AXL', 'TIMD4', 'SLC40A1', 'HMOX1',
      'ANXA1', 'HMGB2', 'NLRP3', 'AIF1', 'SIGLEC1', 'RNASE2'
    )
  )

  save_dotplot(
    sce_obj,
    features = groups$interested,
    assay = assay,
    filename = paste0(prefix, pro, 'DotPlot_check_interested_markers_by_clusters_MonoMacro_cells.pdf'),
    width = 16,
    height = 24
  )

  save_dotplot(
    sce_obj,
    features = groups$univ_mhc,
    assay = assay,
    filename = paste0(prefix, pro, 'DotPlot_univ_mhc_by_clusters.pdf'),
    width = 9,
    height = 10,
    title = 'Universal makers and MHC gene markers'
  )

  save_dotplot(
    sce_obj,
    features = groups$interferon,
    assay = assay,
    filename = paste0(prefix, pro, 'DotPlot_interferon_markers_by_clusters.pdf'),
    width = 9,
    height = 10,
    title = 'Interferon genes'
  )

  save_dotplot(
    sce_obj,
    features = groups$receptors,
    assay = assay,
    filename = paste0(prefix, pro, "DotPlot_fcReceptor_markers_by_clusters.pdf"),
    width = 9,
    height = 10,
    title = "Fc Receptors, Scavenger Receptors, Don't eat me and MMP genes"
  )

  save_dotplot(
    sce_obj,
    features = groups$s100_cytokine,
    assay = assay,
    filename = 'DotPlot_s100_cytokine_markers_by_clusters.pdf',
    width = 9,
    height = 10,
    title = 'S100 gene family, cytokines and other markers'
  )
}

if (length(target_indices) > 0) {
  mono_idx <- target_indices[1]
  with_dir(fname[mono_idx], {
    sce_mono <- readRDS(paste0(fname[mono_idx], ".sce.sub.harmony.Rds"))
    plot_monocyte_markers(sce_mono, prefix = fname[mono_idx], assay = DefaultAssay(sce_mono))
  })
}

### For B and PlasmaCells
with_dir("Plasmacells", {
  sce <- readRDS("Plasmacellssubcluster_best_resolution.Rds")

  genes_to_check <- c(
    'MS4A1', 'CD79B', 'CD79A', 'CD22', 'IGHM', 'IGHD', 'IGLC3', 'CD19', 'CD38', 'SDC1', 'TNFRSF17',
    'SLAMF7', 'GPRC5D', 'CD24', 'GAS6', 'GAS7', 'IGF1', 'CD27'
  )

  save_dotplot(
    sce,
    features = genes_to_check,
    assay = "SCT",
    filename = 'DotPlot_check_interested_markers_for_B_Plasma_cells.pdf',
    width = 10,
    height = 8
  )
})
