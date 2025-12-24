# Step 14: Compare CellChat-derived cell-cell communication across multiple groups

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
  dataset_column = Sys.getenv("CELLCHAT_DATASET_COLUMN", "group"),
  celltype_column = Sys.getenv("CELLCHAT_CELLTYPE_COLUMN", "maincelltype"),
  species = tolower(Sys.getenv("CELLCHAT_SPECIES", "human")),
  min_cells = as.numeric(Sys.getenv("CELLCHAT_MIN_CELLS", "10")),
  do_parallel = tolower(Sys.getenv("CELLCHAT_PARALLEL", "true")) == "true",
  workers = as.numeric(Sys.getenv("CELLCHAT_WORKERS", max(1, future::availableCores() - 1))),
  top_pathways = as.numeric(Sys.getenv("CELLCHAT_TOP_PATHWAYS", "20"))
)

modified_dir <- file.path(config$base_dir, "4-modifiedcelltype")
cellchat_dir <- file.path(config$base_dir, "5-cellchat", "comparisons")

for (dir_path in c(modified_dir, cellchat_dir)) {
  dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)
}

input_rds <- file.path(modified_dir, config$seurat_file)
if (!file.exists(input_rds)) {
  stop("Input Seurat object not found: ", input_rds)
}

sce <- readRDS(input_rds)
meta <- sce@meta.data

if (!config$dataset_column %in% colnames(meta)) {
  stop("Dataset grouping column not found in metadata: ", config$dataset_column)
}
if (!config$celltype_column %in% colnames(meta)) {
  stop("Cell-type grouping column not found in metadata: ", config$celltype_column)
}

groups <- na.omit(unique(meta[[config$dataset_column]]))
if (length(groups) < 2) {
  stop("Need at least two groups to compare. Found: ", paste(groups, collapse = ", "))
}

if (config$do_parallel) {
  future::plan("multisession", workers = config$workers)
  on.exit(future::plan("sequential"), add = TRUE)
}

initialize_cellchat <- function(seurat_obj, celltype_column, species, min_cells) {
  cellchat_obj <- createCellChat(object = seurat_obj, group.by = celltype_column)

  species <- match.arg(species, choices = c("human", "mouse"))
  if (species == "human") {
    cellchat_obj@DB <- CellChatDB.human
    cellchat_obj <- projectData(cellchat_obj, PPI.human)
  } else {
    cellchat_obj@DB <- CellChatDB.mouse
  }

  cellchat_obj <- subsetData(cellchat_obj)
  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
  cellchat_obj <- computeCommunProb(cellchat_obj)
  cellchat_obj <- filterCommunication(cellchat_obj, min.cells = min_cells)
  cellchat_obj <- computeCommunProbPathway(cellchat_obj)
  cellchat_obj <- aggregateNet(cellchat_obj)
  cellchat_obj <- netAnalysis_computeCentrality(cellchat_obj, slot.name = "netP")

  return(cellchat_obj)
}

save_pdf <- function(file, expr, width = 8, height = 8) {
  pdf(file, width = width, height = height)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

cellchat_list <- list()
for (grp in groups) {
  grp_cells <- rownames(meta[meta[[config$dataset_column]] == grp, , drop = FALSE])
  grp_seurat <- subset(sce, cells = grp_cells)
  message("Processing group: ", grp, " (", length(grp_cells), " cells)")
  cellchat_list[[as.character(grp)]] <- initialize_cellchat(
    seurat_obj = grp_seurat,
    celltype_column = config$celltype_column,
    species = config$species,
    min_cells = config$min_cells
  )
}

saveRDS(cellchat_list, file.path(cellchat_dir, paste0("cellchat_list_by_", config$dataset_column, ".rds")))

communication_tables <- bind_rows(lapply(names(cellchat_list), function(grp) {
  df <- subsetCommunication(cellchat_list[[grp]])
  df$dataset <- grp
  df
}))
write.csv(
  communication_tables,
  file = file.path(cellchat_dir, paste0("cellchat_communication_by_", config$dataset_column, ".csv")),
  row.names = FALSE
)

cellchat_combined <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
saveRDS(cellchat_combined, file.path(cellchat_dir, paste0("cellchat_combined_", config$dataset_column, ".rds")))

bubble_file <- file.path(cellchat_dir, paste0("cellchat_bubble_comparison_", config$dataset_column, ".pdf"))
save_pdf(
  bubble_file,
  netVisual_bubble(
    cellchat_combined,
    sources.use = NULL,
    targets.use = NULL,
    signaling = head(rankNet(cellchat_combined, mode = "comparison")$pathways, config$top_pathways),
    comparison = seq_along(cellchat_list),
    angle.x = 45,
    remove.isolate = TRUE
  )
)

heatmap_count_file <- file.path(cellchat_dir, paste0("cellchat_heatmap_count_comparison_", config$dataset_column, ".pdf"))
save_pdf(
  heatmap_count_file,
  netVisual_heatmap(cellchat_combined, measure = "count", color.heatmap = "Blues")
)

heatmap_strength_file <- file.path(cellchat_dir, paste0("cellchat_heatmap_strength_comparison_", config$dataset_column, ".pdf"))
save_pdf(
  heatmap_strength_file,
  netVisual_heatmap(cellchat_combined, measure = "weight", color.heatmap = "Reds")
)

centrality_file <- file.path(cellchat_dir, paste0("cellchat_centrality_comparison_", config$dataset_column, ".pdf"))
save_pdf(
  centrality_file,
  netAnalysis_signalingRole_scatter(cellchat_combined)
)

for (i in seq_along(groups)) {
  for (j in seq_along(groups)) {
    if (i < j) {
      count_plot <- compareInteractions(cellchat_combined, show.legend = FALSE, group = c(i, j))
      ggsave(
        filename = file.path(cellchat_dir, paste0("interaction_counts_", groups[i], "_vs_", groups[j], ".pdf")),
        plot = count_plot,
        width = 6,
        height = 5
      )

      strength_plot <- compareInteractions(cellchat_combined, show.legend = FALSE, group = c(i, j), measure = "weight")
      ggsave(
        filename = file.path(cellchat_dir, paste0("interaction_strength_", groups[i], "_vs_", groups[j], ".pdf")),
        plot = strength_plot,
        width = 6,
        height = 5
      )
    }
  }
}

summary_csv <- file.path(cellchat_dir, paste0("cellchat_summary_counts_strength_", config$dataset_column, ".csv"))
summary_df <- data.frame(
  dataset = names(cellchat_list),
  interaction_count = vapply(cellchat_list, function(obj) sum(obj@net$count), numeric(1)),
  interaction_strength = vapply(cellchat_list, function(obj) sum(obj@net$weight), numeric(1))
)
write.csv(summary_df, summary_csv, row.names = FALSE)

message("Multi-group CellChat comparison completed. Outputs saved to: ", cellchat_dir)
