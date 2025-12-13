# Step 9: map corrected subcell and main cell types back to the full object

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

base_dir <- "/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results"
celltype_dir <- file.path(base_dir, "3-celltype")
modified_dir <- file.path(base_dir, "4-modifiedcelltype")

load_object <- function(path) {
  if (!file.exists(path)) stop(path, " not found")
  readRDS(path)
}

prioritize_columns <- function(df, columns) {
  values <- lapply(columns, function(col) {
    v <- df[[col]]
    v[v == ""] <- NA_character_
    v
  })
  Reduce(function(x, y) ifelse(!is.na(x), x, y), values)
}

sce.epcells <- load_object(file.path(celltype_dir, "Epithelialcells", "Epithelialcells.Rds"))
sce.epcells$subcelltype <- "Epithelial cells"
sce.epcells$maincelltype <- "Epithelial cells"

sce.fibro <- load_object(file.path(celltype_dir, "Fibroblasts", "Fibroblasts.Rds"))
sce.fibro$subcelltype <- "Fibroblasts"
sce.fibro$maincelltype <- "Fibroblasts"

sce.gran <- load_object(file.path(celltype_dir, "Neutrophils", "Neutrophils.Rds"))
sce.gran$subcelltype <- "Granulocytes"
sce.gran$maincelltype <- "Granulocytes"

sce.mast <- load_object(file.path(celltype_dir, "Mastcells", "sce.Mastcells_w_correction.Rds"))

sce.b <- load_object(file.path(modified_dir, "sce.Bcells_w_correction.Rds"))
sce.t <- load_object(file.path(modified_dir, "sce.T_w_correction.Rds"))
sce.dc <- load_object(file.path(modified_dir, "sce.DCs_w_correction.Rds"))
sce.monomacro <- load_object(file.path(modified_dir, "sce.MonoMacro_w_correction.Rds"))

sce.all <- load_object(file.path(celltype_dir, "sce.all_w_cell_type_anno.Rds"))
sce.all$group <- factor(sce.all$group, levels = c("PBS", "DD_mGE", "DD_mIL12", "DD_mGE12", "TD_mGE12"))

sce.all@meta.data$bsubcelltype <- sce.b@meta.data$newsubcelltype[match(rownames(sce.all@meta.data), rownames(sce.b@meta.data))]
sce.all@meta.data$tsubcelltype <- sce.t@meta.data$newsubcelltype[match(rownames(sce.all@meta.data), rownames(sce.t@meta.data))]
sce.all@meta.data$dcsubcelltype <- sce.dc@meta.data$newsubcelltype[match(rownames(sce.all@meta.data), rownames(sce.dc@meta.data))]
sce.all@meta.data$monosubcelltype <- sce.monomacro@meta.data$newsubcelltype[match(rownames(sce.all@meta.data), rownames(sce.monomacro@meta.data))]
sce.all@meta.data$gransubcelltype <- sce.gran@meta.data$subcelltype[match(rownames(sce.all@meta.data), rownames(sce.gran@meta.data))]
sce.all@meta.data$fibrosubcelltype <- sce.fibro@meta.data$subcelltype[match(rownames(sce.all@meta.data), rownames(sce.fibro@meta.data))]
sce.all@meta.data$episubcelltype <- sce.epcells@meta.data$subcelltype[match(rownames(sce.all@meta.data), rownames(sce.epcells@meta.data))]
sce.all@meta.data$mastsubcelltype <- sce.mast@meta.data$subcelltype[match(rownames(sce.all@meta.data), rownames(sce.mast@meta.data))]

subcell_columns <- c(
  "bsubcelltype", "tsubcelltype", "dcsubcelltype", "monosubcelltype",
  "gransubcelltype", "fibrosubcelltype", "episubcelltype", "mastsubcelltype", "subcelltype"
)

sce.all$newsubcelltype <- prioritize_columns(sce.all@meta.data, subcell_columns)

sce.all@meta.data$bmaincelltype <- sce.b@meta.data$maincelltype[match(rownames(sce.all@meta.data), rownames(sce.b@meta.data))]
sce.all@meta.data$tmaincelltype <- sce.t@meta.data$maincelltype[match(rownames(sce.all@meta.data), rownames(sce.t@meta.data))]
sce.all@meta.data$dcmaincelltype <- sce.dc@meta.data$maincelltype[match(rownames(sce.all@meta.data), rownames(sce.dc@meta.data))]
sce.all@meta.data$monomaincelltype <- sce.monomacro@meta.data$maincelltype[match(rownames(sce.all@meta.data), rownames(sce.monomacro@meta.data))]
sce.all@meta.data$granmaincelltype <- sce.gran@meta.data$maincelltype[match(rownames(sce.all@meta.data), rownames(sce.gran@meta.data))]
sce.all@meta.data$fibromaincelltype <- sce.fibro@meta.data$maincelltype[match(rownames(sce.all@meta.data), rownames(sce.fibro@meta.data))]
sce.all@meta.data$epimaincelltype <- sce.epcells@meta.data$maincelltype[match(rownames(sce.all@meta.data), rownames(sce.epcells@meta.data))]
sce.all@meta.data$mastmaincelltype <- sce.mast@meta.data$maincelltype[match(rownames(sce.all@meta.data), rownames(sce.mast@meta.data))]

maincell_columns <- c(
  "bmaincelltype", "tmaincelltype", "dcmaincelltype", "monomaincelltype",
  "granmaincelltype", "fibromaincelltype", "epimaincelltype", "mastmaincelltype", "maincelltype"
)

sce.all$maincelltype <- prioritize_columns(sce.all@meta.data, maincell_columns)

DimPlot(sce.all, raster = FALSE, label = TRUE)
Idents(sce.all) <- sce.all$maincelltype
DimPlot(sce.all, raster = FALSE, label = TRUE)

saveRDS(sce.all, file.path(modified_dir, "sce.all_w_corrected_new_subcelltype_main_celltype.Rds"))

sce.new <- subset(sce.all, subset = maincelltype != "Unknown")
sce.new$subcelltype <- sce.new$newsubcelltype
saveRDS(sce.new, file.path(modified_dir, "sce.new_w_corrected_main_celltype_and_new_subtype.Rds"))

Idents(sce.new) <- sce.new$subcelltype
p <- DimPlot(sce.new, reduction = "umap", label = TRUE, group.by = "subcelltype")
ggsave(
  filename = file.path(modified_dir, "UMAP_subcelltype_A4_cairo.pdf"),
  plot = LabelClusters(plot = p, id = "subcelltype", repel = TRUE),
  device = cairo_pdf, width = 11.69, height = 8.27, dpi = 300
)

saveRDS(sce.new, file.path(modified_dir, "sce.all_w_corrected_new_subcelltype_main_celltype_wo_unknown.Rds"))
