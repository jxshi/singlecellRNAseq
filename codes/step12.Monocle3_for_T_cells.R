### ---------------------------
### Monocle3 pseudotime for Cd4/Cd8 T cells
### ---------------------------
### This script runs monocle3-based trajectories for Cd4 and Cd8 T cells.
### It expects a Seurat object that already contains refined T-cell subtypes
### (e.g., from step8.1_Assign_T_newsubcelltype.R) in the metadata column
### `newsubcelltype`.

rm(list = ls())
options(stringsAsFactors = FALSE)

# -----------------------------
# 0. Setup
# -----------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(dplyr)
})

# Set working directory to where intermediate T-cell objects are stored.
# Adjust the path to your environment before running.
setwd("/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype/newTcells/")

# Load the unified T-cell object with corrected batch effects.
# If you already have split Cd4/Cd8 objects, skip the subsetting section below.
sce.t <- readRDS(
  "/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype/sce.T_w_correction.Rds"
)

# Ensure celltype information exists
sce.t$celltype <- sce.t$newsubcelltype
Idents(sce.t) <- "celltype"

# Subset Cd4 and Cd8 T-cell populations
sce.t.cd4 <- subset(
  sce.t,
  idents = c(
    "Activated Th2", "Cd4 T early activated", "Cd4 Tcm", "Cd4 Tem",
    "Proliferating Treg", "Tfh", "Th1", "Th2", "Tn", "Treg"
  )
)

sce.t.cd8 <- subset(
  sce.t,
  idents = c(
    "Cd8 Teff", "Cd8 Tex", "Cd8 Tm", "Cd8 Tpex", "IFN-responsive Cd8 Teff",
    "Ly6c+ Cd8 Teff", "Proliferating Cd8 T", "Proliferating Cd8 Teff", "Tn"
  )
)

saveRDS(sce.t.cd4, "sce.t.Cd4.Rds")
saveRDS(sce.t.cd8, "sce.t.Cd8.Rds")

# -----------------------------
# 1. Helper: run monocle3 pseudotime with Tn root cluster
# -----------------------------

run_monocle3_pseudotime <- function(seu_subset, root_cluster_label = "Tn", num_dim = 50) {
  stopifnot(is(seu_subset, "Seurat"))

  cds <- as.cell_data_set(seu_subset)

  # Carry over celltype annotation
  colData(cds)$celltype <- seu_subset$celltype

  # Preprocess and reduce dimensionality
  cds <- preprocess_cds(cds, num_dim = num_dim)

  # Optional batch alignment if a batch column exists
  if ("batch" %in% colnames(colData(cds))) {
    cds <- align_cds(cds, alignment_group = "batch")
  }

  cds <- reduce_dimension(cds, reduction_method = "UMAP")
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cds <- learn_graph(cds, use_partition = FALSE)

  # Identify root cells by cluster label (e.g., naive T cells marked as "Tn")
  cluster_info <- clusters(cds)
  if (!root_cluster_label %in% cluster_info) {
    # If clusters are encoded as numeric, map levels to labels in metadata
    if (root_cluster_label %in% cds@colData$celltype) {
      root_cells <- colnames(cds)[cds@colData$celltype == root_cluster_label]
    } else {
      stop("root_cluster_label not found in clusters or celltype metadata: ", root_cluster_label)
    }
  } else {
    root_cells <- colnames(cds)[cluster_info == root_cluster_label]
  }

  cds <- order_cells(cds, root_cells = root_cells)
  cds
}

# -----------------------------
# 2. Run trajectories for Cd4 and Cd8 T cells
# -----------------------------

cds_cd4 <- run_monocle3_pseudotime(sce.t.cd4, root_cluster_label = "Tn", num_dim = 50)
cds_cd8 <- run_monocle3_pseudotime(sce.t.cd8, root_cluster_label = "Tn", num_dim = 50)

saveRDS(cds_cd4, "cds_cd4_monocle3.rds")
saveRDS(cds_cd8, "cds_cd8_monocle3.rds")

# -----------------------------
# 3. Visualization and differential tests
# -----------------------------

pdf("monocle3_cd4_cd8_pseudotime_plots.pdf", width = 12, height = 6)

par(mfrow = c(2, 2))
plot_cells(cds_cd4, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
plot_cells(cds_cd4, color_cells_by = "celltype", show_trajectory_graph = TRUE)
plot_cells(cds_cd8, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
plot_cells(cds_cd8, color_cells_by = "celltype", show_trajectory_graph = TRUE)

dev.off()

# Differential expression along pseudotime
# Note: graph_test can be slow; adjust cores to your environment.
deg_cd4 <- graph_test(cds_cd4, neighbor_graph = "principal_graph", cores = 4)
deg_cd8 <- graph_test(cds_cd8, neighbor_graph = "principal_graph", cores = 4)

write.csv(deg_cd4, "deg_cd4_pseudotime.csv")
write.csv(deg_cd8, "deg_cd8_pseudotime.csv")

message("Monocle3 pseudotime completed for Cd4 and Cd8 T cells with Tn as the root cluster.")
