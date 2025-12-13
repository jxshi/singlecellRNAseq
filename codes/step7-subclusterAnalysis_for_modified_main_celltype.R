# Step 7: split objects by corrected main cell type and save per-celltype files

load_required_packages <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
  }
  invisible(lapply(pkgs, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
}

load_required_packages(c("Seurat", "ggplot2", "clustree", "cowplot", "dplyr"))

base_dir <- "/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype"
setwd(base_dir)

assert_file <- function(path) {
  if (!file.exists(path)) stop(path, " not found")
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

save_subset <- function(object, celltype_name, folder_name) {
  ensure_dir(folder_name)
  subset_obj <- subset(object, idents = celltype_name)
  saveRDS(subset_obj, file = file.path(folder_name, paste0(folder_name, ".Rds")))
  write.table(
    table(subset_obj$celltype, subset_obj$group),
    file = file.path(folder_name, paste0(folder_name, "_celltype_by_group.tsv")),
    sep = "\t", quote = FALSE, col.names = NA
  )
  subset_obj
}

sce_path <- "sce.new_w_corrected_main_celltype.Rds"
assert_file(sce_path)
sce_all <- readRDS(sce_path)

sce_all.list <- SplitObject(sce_all, split.by = "maincelltype")
celltype_names <- names(sce_all.list)
folder_names <- gsub(" ", "", celltype_names)

saveRDS(celltype_names, "celltype_name.Rds")
saveRDS(folder_names, "celltype_file_folder.Rds")

invisible(mapply(save_subset, celltype_names, folder_names, MoreArgs = list(object = sce_all)))

# Quick sanity check on macrophage markers before downstream work
qc_markers <- c("Cd14", "Cd16", "Marco")
qc_plot <- FeaturePlot(sce_all, qc_markers)
ggsave("featureplot_modified_main_celltype_qc.pdf", qc_plot, width = 10, height = 8)
