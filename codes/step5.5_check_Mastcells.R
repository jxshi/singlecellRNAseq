# Step 5.5: check Mast cell subclusters

load_required_packages <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
  }
  invisible(lapply(pkgs, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
}

load_required_packages(c("Seurat", "dplyr"))

working_dir <- "/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/3-celltype/Mastcells"
setwd(working_dir)

assert_file <- function(path) {
  if (!file.exists(path)) {
    stop(path, " not found")
  }
}

apply_cluster_labels <- function(obj, mapping, column) {
  obj[[column]] <- "NA"
  for (i in seq_len(nrow(mapping))) {
    obj@meta.data[obj$seurat_clusters == mapping$ClusterID[i], column] <- mapping[[column]][i]
  }
  obj
}

write_proportion_tables <- function(obj, column, prefix = column) {
  clty_per_ident <- table(obj@meta.data[[column]], obj@meta.data$orig.ident)
  write.table(
    clty_per_ident,
    paste0(prefix, "_per_orig.ident.tsv"),
    quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA
  )

  clty_per_group <- table(obj@meta.data[[column]], obj@meta.data$group)
  write.table(
    clty_per_group,
    paste0(prefix, "_per_group.tsv"),
    quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA
  )
}

assert_file("Mastcells.sce.sub.harmony.Rds")
sce.sub <- readRDS("Mastcells.sce.sub.harmony.Rds")

subcelltype <- data.frame(ClusterID = 0:9, subcelltype = "na")
subcelltype[subcelltype$ClusterID %in% c(0), 2] <- "Unknown"
subcelltype[subcelltype$ClusterID %in% c(1), 2] <- "CTMC"
subcelltype[subcelltype$ClusterID %in% c(2), 2] <- "MMC"
subcelltype[subcelltype$ClusterID %in% c(3), 2] <- "Unknown"
subcelltype[subcelltype$ClusterID %in% c(4), 2] <- "Unknown"
subcelltype[subcelltype$ClusterID %in% c(5), 2] <- "Unknown"
subcelltype[subcelltype$ClusterID %in% c(6), 2] <- "CTMC"
subcelltype[subcelltype$ClusterID %in% c(7), 2] <- "CTMC"
subcelltype[subcelltype$ClusterID %in% c(8), 2] <- "Unknown"
subcelltype[subcelltype$ClusterID %in% c(9), 2] <- "Unknown"

sce.sub <- apply_cluster_labels(sce.sub, subcelltype, "subcelltype")
write_proportion_tables(sce.sub, "subcelltype")
Idents(sce.sub) <- sce.sub$subcelltype

maincelltype <- data.frame(ClusterID = 0:9, maincelltype = "na")
maincelltype[maincelltype$ClusterID %in% c(0), 2] <- "Unknown"
maincelltype[maincelltype$ClusterID %in% c(1), 2] <- "Mast cells"
maincelltype[maincelltype$ClusterID %in% c(2), 2] <- "Mast cells"
maincelltype[maincelltype$ClusterID %in% c(3), 2] <- "Unknown"
maincelltype[maincelltype$ClusterID %in% c(4), 2] <- "Unknown"
maincelltype[maincelltype$ClusterID %in% c(5), 2] <- "Unknown"
maincelltype[maincelltype$ClusterID %in% c(6), 2] <- "Mast cells"
maincelltype[maincelltype$ClusterID %in% c(7), 2] <- "Mast cells"
maincelltype[maincelltype$ClusterID %in% c(8), 2] <- "Unknown"
maincelltype[maincelltype$ClusterID %in% c(9), 2] <- "Unknown"

sce.sub <- apply_cluster_labels(sce.sub, maincelltype, "maincelltype")
write_proportion_tables(sce.sub, "maincelltype")
Idents(sce.sub) <- sce.sub$maincelltype

sce.mast <- sce.sub

table(sce.mast$group, sce.mast$maincelltype)

saveRDS(sce.mast, "sce.Mastcells_w_correction.Rds")
