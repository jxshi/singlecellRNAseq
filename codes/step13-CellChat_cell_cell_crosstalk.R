# Step 13: CellChat-based cell-cell communication analysis

load_required_packages <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
  }
  invisible(lapply(pkgs, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
}

load_required_packages(c("Seurat", "CellChat", "patchwork", "ggplot2", "dplyr"))

base_dir <- "/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results"
modified_dir <- file.path(base_dir, "4-modifiedcelltype")
cellchat_dir <- file.path(base_dir, "5-cellchat")

dir.create(cellchat_dir, showWarnings = FALSE, recursive = TRUE)

input_rds <- file.path(modified_dir, "sce.all_w_corrected_new_subcelltype_main_celltype_wo_unknown.Rds")
if (!file.exists(input_rds)) {
  stop("Input Seurat object not found: ", input_rds)
}

sce <- readRDS(input_rds)

# Choose the grouping column for communication analysis
# Options: "maincelltype" or "subcelltype"
grouping_column <- "maincelltype"

if (!grouping_column %in% colnames(sce@meta.data)) {
  stop("Grouping column not found in metadata: ", grouping_column)
}

# Select species database ("human" or "mouse")
species <- "human"

cellchat <- createCellChat(object = sce, group.by = grouping_column)

if (species == "human") {
  cellchat@DB <- CellChatDB.human
  cellchat <- projectData(cellchat, PPI.human)
} else if (species == "mouse") {
  cellchat@DB <- CellChatDB.mouse
} else {
  stop("Unsupported species: ", species)
}

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file.path(cellchat_dir, "cellchat_maincelltype.rds"))

# Save the inferred communication table
communication_table <- subsetCommunication(cellchat)
write.csv(
  communication_table,
  file = file.path(cellchat_dir, "cellchat_communication_table.csv"),
  row.names = FALSE
)

pdf(file.path(cellchat_dir, "cellchat_net_circle_interactions.pdf"), width = 8, height = 8)
netVisual_circle(
  cellchat@net$count,
  vertex.weight = as.numeric(table(cellchat@idents)),
  weight.scale = TRUE,
  label.edge = FALSE
)
dev.off()

pdf(file.path(cellchat_dir, "cellchat_net_circle_strength.pdf"), width = 8, height = 8)
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = as.numeric(table(cellchat@idents)),
  weight.scale = TRUE,
  label.edge = FALSE
)
dev.off()

pdf(file.path(cellchat_dir, "cellchat_net_pathways.pdf"), width = 10, height = 8)
netAnalysis_signalingRole_network(cellchat, signaling = "all")
dev.off()
