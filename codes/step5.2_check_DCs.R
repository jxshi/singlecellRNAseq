# Step 5.2: check DC subclusters

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
  "Seurat", "tidyverse", "ggplot2", "sctransform", "glmGamPoi", "DoubletFinder", "patchwork",
  "clusterProfiler", "org.Hs.eg.db", "harmony", "celldex", "RColorBrewer", "future", "parallel"
))

source("~/software/functions/custom_seurat_functions.R")
source("~/software/functions/PropPlot.R")
source("~/software/functions/SubClusterPropPlot.R")
source("~/software/functions/convertHumanGeneList.R")

options(future.globals.maxSize = 891289600)
options(future.seed = TRUE)

dotplot_theme <- theme(
  axis.text.x = element_text(angle = 45, hjust = 1.2, vjust = 1.1, size = 8),
  strip.text  = element_text(size = 8)
)

save_dotplot <- function(object, genes, filename, assay = "RNA", width = 14, height = 10, title = NULL) {
  plot <- DotPlot(object, features = genes, assay = assay) + dotplot_theme
  if (!is.null(title)) plot <- plot + ggtitle(title)
  ggsave(filename, plot = plot, width = width, height = height)
  plot
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

working_dir <- "/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/3-celltype/DCs"
if (dir.exists(working_dir)) setwd(working_dir)

prefix <- "DCs"
dataset_path <- "DCs.sce.sub.harmony.Rds"
if (!file.exists(dataset_path)) stop(dataset_path, " not found")

sce <- readRDS(dataset_path)

# DC subtype marker genes (mouse; from your Excel)
genes_dc <- list(
  DC_cDC1   = c("Xcr1", "Wdfy4", "Clec9a"),
  DC_cDC2   = c("Cd209a", "Clec10a", "Mgl2"),
  `DC_成熟DC` = c("Fscn1", "Ccl22", "Nudt11"),
  DC_pDC    = c("Irt7", "Spib", "Bst2", "Siglech", "Ccr9", "CD209a", "Upb1", "Rnd3", "Runx2")
)

p_all_markers <- DotPlot(sce, features = genes_dc, assay = "RNA") + dotplot_theme

p_all_markers
save_dotplot(sce, genes_dc, "check_all_cell_markers_for_DCs.pdf", width = 14, height = 10)

subcelltype <- data.frame(
  ClusterID = 0:12,
  subcelltype = c(
    "cDC2", "moDC", "cDC1", "pDC", "mregDC", "moDC", "cDC2", "Unknown", "cDC1", "cDC2", "mregDC", "pDC", "pDC"
  )
)

sce <- apply_cluster_labels(sce, subcelltype, "subcelltype")
write_proportion_tables(sce, "subcelltype", prefix)
Idents(sce) <- sce$subcelltype

maincelltype <- data.frame(
  ClusterID = 0:12,
  maincelltype = c(
    "DCs", "DCs", "DCs", "DCs", "DCs", "DCs", "DCs", "Unknown", "DCs", "DCs", "DCs", "DCs", "DCs"
  )
)

sce <- apply_cluster_labels(sce, maincelltype, "maincelltype")
write_proportion_tables(sce, "maincelltype", prefix)
Idents(sce) <- sce$maincelltype

saveRDS(sce, "sce.DCs_w_correction.Rds")

# Proportion tables for quick review
table(sce$group, sce$maincelltype)
table(sce$group, sce$subcelltype)
