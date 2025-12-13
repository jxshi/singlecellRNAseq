#################### Step 3.1 Cell Type Recognition ####################
if (!dir.exists("../3-celltype")) {
  dir.create("../3-celltype")
}
setwd("../3-celltype")

silent_load <- function(pkgs) {
  invisible(lapply(pkgs, function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
      stop("Package not installed: ", pkg)
    }
  }))
}

run_balloonplot <- function(tab, file, width = 16, height = 10) {
  pdf(file, width = width, height = height)
  balloonplot(tab)
  dev.off()
}

save_plot <- function(plot_obj, filename, width = 14, height = 10, device = cairo_pdf) {
  ggsave(plot = plot_obj, filename = filename, device = device, width = width, height = height)
}

annotate_with_singleR <- function(sce, ref_expr, ref_main, ref_fine, prefix) {
  main_res <- SingleR(test = sce, assay.type.test = 1, ref = ref_expr, labels = ref_main)
  fine_res <- SingleR(test = sce, assay.type.test = 1, ref = ref_expr, labels = ref_fine)

  message("Completed SingleR annotation for ", prefix)
  list(main = main_res, fine = fine_res)
}

apply_annotations <- function(seurat_obj, res, prefix, cluster_col = "seurat_clusters") {
  meta_main <- paste0(prefix, ".main")
  meta_fine <- paste0(prefix, ".fine")

  seurat_obj@meta.data[[meta_main]] <- res$main$pruned.labels
  seurat_obj@meta.data[[meta_fine]] <- res$fine$pruned.labels

  main_tab <- table(seurat_obj@meta.data[[meta_main]], seurat_obj@meta.data[[cluster_col]])
  fine_tab <- table(seurat_obj@meta.data[[meta_fine]], seurat_obj@meta.data[[cluster_col]])

  list(main_tab = main_tab, fine_tab = fine_tab, seurat_obj = seurat_obj)
}

export_annotation_tables <- function(prefix, main_tab, fine_tab) {
  write.table(main_tab, file = paste0(prefix, ".main.table_for_sce.all.tsv"), sep = "\t", quote = FALSE)
  write.table(fine_tab, file = paste0(prefix, ".fine.table_for_sce.all.tsv"), sep = "\t", quote = FALSE)

  run_balloonplot(main_tab, paste0(prefix, ".main.label_balloonplot.pdf"))
  run_balloonplot(fine_tab, paste0(prefix, ".fine.label_balloonplot.pdf"))
}

silent_load(c("celldex", "SingleR", "Seurat", "gplots", "tidyverse", "patchwork"))

sce.all <- readRDS("../1-QC/sce.all_best_resolution.Rds")

# Use SingleR to annotate cell types to help identify clusters
run_singleR_annotations <- TRUE
if (run_singleR_annotations) {
  sce <- as.SingleCellExperiment(DietSeurat(sce.all))

  references <- c("monaco", "hpca", "dice")

  if ("monaco" %in% references) {
    monaco.ref <- celldex::MonacoImmuneData()
    monaco_res <- annotate_with_singleR(sce, monaco.ref, monaco.ref$label.main, monaco.ref$label.fine, "monaco")
    monaco_tabs <- apply_annotations(sce.all, monaco_res, prefix = "monaco")
    sce.all <- monaco_tabs$seurat_obj
    export_annotation_tables("monaco", monaco_tabs$main_tab, monaco_tabs$fine_tab)
  }

  if ("hpca" %in% references) {
    hpca.ref <- celldex::HumanPrimaryCellAtlasData()
    hpca_res <- annotate_with_singleR(sce, hpca.ref, hpca.ref$label.main, hpca.ref$label.fine, "hpca")
    hpca_tabs <- apply_annotations(sce.all, hpca_res, prefix = "hpca")
    sce.all <- hpca_tabs$seurat_obj
    export_annotation_tables("hpca", hpca_tabs$main_tab, hpca_tabs$fine_tab)
  }

  if ("dice" %in% references) {
    dice.ref <- celldex::DatabaseImmuneCellExpressionData()
    dice_res <- annotate_with_singleR(sce, dice.ref, dice.ref$label.main, dice.ref$label.fine, "dice")
    dice_tabs <- apply_annotations(sce.all, dice_res, prefix = "dice")
    sce.all <- dice_tabs$seurat_obj
    export_annotation_tables("dice", dice_tabs$main_tab, dice_tabs$fine_tab)
  }
}

#################### Step 3.2 Assign cell type to each cluster ####################
# **************** IMPORTANT ****************
# This part has to be done manually.
# Option 1:
run_manual_assignment <- TRUE
if (run_manual_assignment) {
  message("Active cluster counts")
  print(table(sce.all@active.ident))

  celltype <- data.frame(
    ClusterID = 0:28,
    celltype = "na",
    stringsAsFactors = FALSE
  )

  replacements <- list(
    `0` = "T cells",
    `1` = "T cells",
    `2` = "T cells",
    `3` = "T cells",
    `4` = "Mono Macro",
    `5` = "Mono Macro",
    `6` = "Mono Macro",
    `7` = "Mono Macro",
    `8` = "Mono Macro",
    `9` = "B cells",
    `10` = "T cells",
    `11` = "T cells",
    `12` = "T cells",
    `13` = "DCs",
    `14` = "Mono Macro",
    `15` = "T cells",
    `16` = "T cells",
    `17` = "Mono Macro",
    `18` = "Mono Macro",
    `19` = "Epithelial cells",
    `20` = "T cells",
    `21` = "T cells",
    `22` = "Mast cells",
    `23` = "T cells",
    `24` = "DCs",
    `25` = "DCs",
    `26` = "Mono Macro",
    `27` = "Fibroblasts",
    `28` = "Neutrophils"
  )

  for (cluster_id in names(replacements)) {
    celltype[celltype$ClusterID == as.integer(cluster_id), "celltype"] <- replacements[[cluster_id]]
  }

  print(head(celltype))
  print(table(celltype$celltype))

  sce.all@meta.data$celltype <- "NA"
  for (i in seq_len(nrow(celltype))) {
    idx <- sce.all@meta.data$seurat_clusters == celltype$ClusterID[i]
    sce.all@meta.data[idx, "celltype"] <- celltype$celltype[i]
  }

  message("Assigned celltype counts")
  print(table(sce.all@meta.data$celltype))

  cltyPerIdent <- table(sce.all@meta.data$celltype, sce.all@meta.data$orig.ident)
  write.table(cltyPerIdent, "celltype_per_orig.ident.tsv", quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

  cltyPerGroup <- table(sce.all@meta.data$celltype, sce.all@meta.data$group)
  write.table(cltyPerGroup, "celltype_per_group.tsv", quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

  Idents(sce.all) <- sce.all$celltype
}

# Plot proportion plot of cell types
run_proportion_plots <- TRUE
if (run_proportion_plots) {
  source("~/software/functions/PropPlot.R")

  rotate_x <- function(p) {
    p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }

  # Plot proportions of different celltypes per group
  p_group <- PropPlot(sce.all, "group")
  ggsave(rotate_x(p_group), filename = "group_celltype_proportion_plot.pdf", device = cairo_pdf, width = 6, height = 10)

  # Plot proportions of different celltypes per batch
  p_batch <- PropPlot(sce.all, "batch")
  ggsave(rotate_x(p_batch), filename = "batch_celltype_proportion_plot.pdf", device = cairo_pdf, width = 6, height = 10)
}

# Option 2:
run_quick_idents <- FALSE
if (run_quick_idents) {
  sce.all <- RenameIdents(
    sce.all,
    `0` = "Mono/Macro",
    `1` = "B cells",
    `2` = "CD4+ Memory T",
    `3` = "CD8+ NK cells",
    `4` = "NK cells",
    `5` = "Fibroblasts",
    `6` = "Endothelial cells",
    `7` = "Mono/Macro cells",
    `8` = "Fibroblasts",
    `9` = "Myeloid",
    `10` = "Mono/Macro cells",
    `11` = "Myeloid",
    `12` = "Naive B cells",
    `13` = "Plasma B cells",
    `14` = "DC",
    `15` = "NK cells"
  )
  DimPlot(sce.all, label = TRUE)
  saveRDS(sce.all, "sce.all_w_cell_type_anno.Rds")
}
# **************** IMPORTANT ****************
#################### Step 3.3 Check known markers for assigned types of cells ####################
run_marker_review <- TRUE
if (run_marker_review) {
  marker_panel <- c(
    'PTPRC', 'CD3D', 'CD3E', 'CD4', 'CD8A', 'CD19', 'CD79A', 'MS4A1',
    'IGHG1', 'MZB1', 'SDC1',
    'CD68', 'CD163', 'CD14', 'CD206',
    'TPSAB1', 'TPSB2',
    'RCVRN', 'FPR1', 'ITGAM',
    'C1QA', 'C1QB',
    'S100A9', 'S100A8', 'MMP19',
    'LAMP3', 'IDO1', 'IDO2', 'ITGAX', 'CD74', 'FLT3', 'CLEC9A', 'CLEC10A', 'TCF4',
    'CD1E', 'CD1C',
    'KLRB1', 'NCR1',
    'FGF7', 'MME', 'ACTA2',
    'PECAM1', 'VWF',
    'MKI67', 'TOP2A',
    'EPCAM', 'KRT19', 'PROM1', 'ALDH1A1'
  )
  marker_panel <- convertHumanGeneList(marker_panel)

  dot_clusters <- DotPlot(sce.all, features = unique(marker_panel), group.by = "seurat_clusters", assay = 'RNA') +
    coord_flip()
  save_plot(dot_clusters, "check_all_marker_by_seurat_clusters.pdf")

  umap_celltype <- DimPlot(sce.all, reduction = "umap", group.by = "celltype", label = TRUE, raster = FALSE)
  save_plot(umap_celltype, "umap_by_celltype.pdf")

  cluster_tab <- table(sce.all@meta.data$celltype, sce.all@meta.data$seurat_clusters)
  run_balloonplot(cluster_tab, "cluster_celltype_vs_seurat_clusters.pdf")

  rotated_theme <- theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
  dot_celltype <- DotPlot(sce.all, features = unique(marker_panel), assay = 'RNA', group.by = 'celltype') +
    coord_flip() + rotated_theme
  marker_umap <- dot_celltype + umap_celltype
  save_plot(marker_umap, "markers_umap_by_celltype.pdf", width = 16, height = 12, device = "pdf")

  umap_split <- DimPlot(sce.all, reduction = "umap", group.by = "celltype", split.by = "orig.ident", label = TRUE, raster = FALSE)
  save_plot(umap_split, "umap_by_celltype_split_orig.ident.pdf", width = 24, height = 12)

  saveRDS(sce.all, "sce.all_w_cell_type_anno.Rds")
}

run_marker_discovery <- TRUE
if (run_marker_discovery) {
  phe <- sce.all@meta.data
  saveRDS(phe, file = 'phe_by_markers.Rds')

  sce <- sce.all
  Idents(sce) <- sce$celltype
  DefaultAssay(sce) <- "RNA"
  table(Idents(sce))
  sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
  DT::datatable(sce.markers)

  pro <- 'celltype_deg'
  write.csv(sce.markers, file = paste0(pro, '_sce.markers.csv'))

  plot_top_markers <- function(marker_df, top_n, prefix, heatmap_width = 14, heatmap_height = 10, dot_width = 14, dot_height = 10) {
    top_markers <- marker_df %>% group_by(cluster) %>% slice_max(avg_log2FC, n = top_n)
    heatmap_plot <- DoHeatmap(sce, top_markers$gene, size = 3)
    save_plot(heatmap_plot, filename = paste0(prefix, '_DoHeatmap_check_top', top_n, '_markers_by_clusters.pdf'), width = heatmap_width, height = heatmap_height, device = "pdf")

    dot_plot <- DotPlot(sce, features = unique(top_markers$gene), assay = 'RNA') + coord_flip()
    save_plot(dot_plot, filename = paste0(prefix, 'DotPlot_check_top', top_n, '_markers_by_clusters.pdf'), width = dot_width, height = dot_height, device = "pdf")
  }

  plot_top_markers(sce.markers, top_n = 10, prefix = pro, heatmap_width = 14, heatmap_height = 10, dot_width = 14, dot_height = 10)
  plot_top_markers(sce.markers, top_n = 6, prefix = pro, heatmap_width = 18, heatmap_height = 12, dot_width = 18, dot_height = 12)
  plot_top_markers(sce.markers, top_n = 3, prefix = pro)

  saveRDS(sce.markers, file = paste0(pro, 'sce.cell.type.markers.Rds'))
}

# Feature plots of interested marker genes
run_feature_plots <- TRUE
if (run_feature_plots) {
  DefaultAssay(sce) <- "RNA"
  genes_to_check <- c('CD3D','CD4','CD8A','GNLY','IL7R','GZMK','FCGR3A','CD14','CD68','CD163','TNFRSF17','CD19',
                      'CD22','MS4A1','CD79A','IGKC','CD1C','TPM2')
  genes_to_check <- convertHumanGeneList(genes_to_check)

  feature_plots <- FeaturePlot(sce, genes_to_check, raster = FALSE)
  save_plot(feature_plots, "major_cluster_feature_plot.pdf", width = 14, height = 14)
}

#################### Step 3.4 Split data according to cell types ####################
run_split_by_celltype <- TRUE
if (run_split_by_celltype) {
  sce.all.list <- SplitObject(sce.all, split.by = "celltype")
  cname <- names(sce.all.list)
  fname <- gsub(' ', '', names(sce.all.list))

  saveRDS(cname, "celltype_name.Rds")
  saveRDS(fname, "celltype_file_folder.Rds")

  lapply(fname, dir.create, showWarnings = FALSE)

  for (i in seq_along(sce.all.list)) {
    sce.sub <- subset(sce.all, idents = cname[i])
    table(sce.sub$celltype, sce.sub$group)
    saveRDS(sce.sub, file = file.path(fname[i], paste0(fname[i], '.Rds')))
  }
}

### DotPlot for selected genes
run_selected_marker_plots <- TRUE
if (run_selected_marker_plots) {
  interested_markers <- c(
    "CD14", "CD68", "LYZ", "FCGR3B", "CXCR2", "S100A8", "S100A9", "IL7R", "IL32", "CD3D", "CD3G", "CD3E", "CD2",
    "IGHM", "MS4A1", "CD79A", "CD79B", "IGHD", "GNLY", "NKG7", "GZMB", "KLRD1", "FCER1A", "CD1C", "MRC1",
    "GATA1", "GATA2"
  )

  plot_interested <- DotPlot(sce.all, features = unique(interested_markers), group.by = "seurat_clusters", assay = 'RNA') + coord_flip()
  save_plot(plot_interested, "check_interested_marker_by_seurat_clusters.pdf")

  plot_by_celltype <- DotPlot(sce.all, features = unique(interested_markers), group.by = "celltype", assay = 'RNA') + coord_flip()
  save_plot(plot_by_celltype, "check_interested_marker_by_assigned_cell_type.pdf", width = 10, height = 6)

  umap_celltype_final <- DimPlot(sce.all, reduction = "umap", group.by = "celltype", label = TRUE)
  save_plot(umap_celltype_final, 'umap_by_celltype.pdf')
}
