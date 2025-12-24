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

load_required_packages(c("Seurat", "CellChat", "patchwork", "ggplot2", "dplyr", "future"))

config <- list(
  base_dir = "/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results",
  seurat_file = "sce.all_w_corrected_new_subcelltype_main_celltype_wo_unknown.Rds",
  grouping_column = Sys.getenv("CELLCHAT_GROUPING", "maincelltype"),
  species = tolower(Sys.getenv("CELLCHAT_SPECIES", "human")),
  min_cells = as.numeric(Sys.getenv("CELLCHAT_MIN_CELLS", "10")),
  do_parallel = tolower(Sys.getenv("CELLCHAT_PARALLEL", "true")) == "true",
  workers = as.numeric(Sys.getenv("CELLCHAT_WORKERS", max(1, future::availableCores() - 1))),
  top_pathways = as.numeric(Sys.getenv("CELLCHAT_TOP_PATHWAYS", "20"))
)

modified_dir <- file.path(config$base_dir, "4-modifiedcelltype")
cellchat_dir <- file.path(config$base_dir, "5-cellchat")

dir.create(cellchat_dir, showWarnings = FALSE, recursive = TRUE)

input_rds <- file.path(modified_dir, config$seurat_file)
if (!file.exists(input_rds)) {
  stop("Input Seurat object not found: ", input_rds)
}

sce <- readRDS(input_rds)

if (!config$grouping_column %in% colnames(sce@meta.data)) {
  stop("Grouping column not found in metadata: ", config$grouping_column)
}

if (config$do_parallel) {
  future::plan("multisession", workers = config$workers)
  on.exit(future::plan("sequential"), add = TRUE)
}

cellchat <- createCellChat(object = sce, group.by = config$grouping_column)

config$species <- match.arg(config$species, choices = c("human", "mouse"))
if (config$species == "human") {
  cellchat@DB <- CellChatDB.human
  cellchat <- projectData(cellchat, PPI.human)
} else {
  cellchat@DB <- CellChatDB.mouse
}

message("Subset data and detect overexpressed signals ...")
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

message("Compute communication probability ...")
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = config$min_cells)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

saveRDS(cellchat, file.path(cellchat_dir, paste0("cellchat_", config$grouping_column, ".rds")))

communication_table <- subsetCommunication(cellchat)
write.csv(
  communication_table,
  file = file.path(cellchat_dir, paste0("cellchat_communication_table_", config$grouping_column, ".csv")),
  row.names = FALSE
)

aggregate_df <- bind_rows(lapply(names(cellchat@net), function(slot) {
  mat <- cellchat@net[[slot]]
  as.data.frame(as.table(mat)) %>%
    rename(source = Var1, target = Var2, value = Freq) %>%
    mutate(type = slot)
}))
write.csv(aggregate_df, file = file.path(cellchat_dir, paste0("cellchat_network_matrices_", config$grouping_column, ".csv")), row.names = FALSE)

save_pdf <- function(file, expr, width = 8, height = 8) {
  pdf(file, width = width, height = height)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

save_pdf(
  file.path(cellchat_dir, paste0("cellchat_group_size_", config$grouping_column, ".pdf")),
  netVisual_circle(cellchat@net$count, vertex.weight = as.numeric(table(cellchat@idents)), weight.scale = TRUE, label.edge = FALSE)
)

save_pdf(
  file.path(cellchat_dir, paste0("cellchat_interaction_strength_", config$grouping_column, ".pdf")),
  netVisual_circle(cellchat@net$weight, vertex.weight = as.numeric(table(cellchat@idents)), weight.scale = TRUE, label.edge = FALSE)
)

save_pdf(
  file.path(cellchat_dir, paste0("cellchat_pathway_network_", config$grouping_column, ".pdf")),
  netAnalysis_signalingRole_network(cellchat, signaling = "all")
)

save_pdf(
  file.path(cellchat_dir, paste0("cellchat_pathway_heatmap_", config$grouping_column, ".pdf")),
  netVisual_heatmap(cellchat, color.heatmap = "Reds")
)

pathway_strength <- sapply(cellchat@netP$prob, function(mat) sum(mat, na.rm = TRUE))
pathway_use <- head(names(sort(pathway_strength, decreasing = TRUE)), config$top_pathways)
save_pdf(
  file.path(cellchat_dir, paste0("cellchat_bubble_top_pathways_", config$grouping_column, ".pdf")),
  netVisual_bubble(cellchat, sources.use = NULL, targets.use = NULL, signaling = pathway_use, comparison = c(1, 1), angle.x = 45, remove.isolate = TRUE)
)

save_pdf(
  file.path(cellchat_dir, paste0("cellchat_centrality_roles_", config$grouping_column, ".pdf")),
  netAnalysis_signalingRole_scatter(cellchat)
)

message("CellChat analysis completed. Outputs saved to: ", cellchat_dir)
