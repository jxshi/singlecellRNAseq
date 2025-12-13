# Step 5.3: check Mono/Macrophage subclusters

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
  "Seurat", "tidyverse", "ggplot2", "sctransform", "glmGamPoi", "DoubletFinder",
  "patchwork", "clusterProfiler", "org.Hs.eg.db", "harmony", "celldex",
  "RColorBrewer", "future", "parallel"
))

source("~/software/functions/custom_seurat_functions.R")
source("~/software/functions/PropPlot.R")
source("~/software/functions/SubClusterPropPlot.R")
source("~/software/functions/convertHumanGeneList.R")

options(future.globals.maxSize = 891289600)
options(future.seed = TRUE)

working_dir <- "/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/3-celltype/MonoMacro"
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

save_dotplot <- function(object, genes, filename, assay = "RNA", width = 14, height = 10) {
  plot <- DotPlot(object, features = genes, assay = assay) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1.2, vjust = 1.1, size = 8),
      strip.text = element_text(size = 8)
    )
  ggsave(filename, plot = plot, width = width, height = height)
  plot
}

assert_file("MonoMacro.sce.sub.harmony.Rds")
sce <- readRDS("MonoMacro.sce.sub.harmony.Rds")

# Curated mouse macrophage marker sets (usable with Seurat::DotPlot etc.)
genes_mac_curated <- list(
  Mono_Ly6Chi        = c("Ly6c2", "Ccr2", "S100a8", "S100a9", "Plac8", "Lcn2"),
  Mono_Ly6Clo        = c("Nr4a1", "Cx3cr1", "Itgal", "Selplg", "Fcgr3"),
  Mac_TRM_core       = c("Adgre1", "Lyz2", "Mafb", "Apoe", "Lgals3", "C1qa", "C1qb", "C1qc", "Mertk"),
  M1_inflammatory    = c("Nos2", "Il1b", "Tnf", "Il12b", "Ccl2", "Cxcl10", "Cd86", "Stat1", "Irf5"),
  M2a_alternative    = c("Mrc1", "Retnla", "Chil3", "Arg1", "Ccl22", "Il10ra"),
  M2c_deactivated    = c("Cd163", "Tgfb1", "Il10"),
  TAM_general        = c("Cd68", "Csf1r", "Maf", "Axl", "Vegfa", "Mmp9"),
  TAM_MHCII_high     = c("H2-Aa", "H2-Ab1", "Cd74", "Ciita", "Irf1"),
  ISG_macrophage     = c("Isg15", "Ifit1", "Ifit3", "Rsad2", "Oasl1", "Irf7"),
  LAM_Trem2_plus     = c("Trem2", "Cd9", "Gpnmb", "Lpl", "Lipa", "Fabp5", "Cd63", "Ctsb", "Ctsd"),
  SPP1_plus_TAM      = c("Spp1", "Fn1"),
  Kupffer            = c("Clec4f", "Vsig4", "Timd4", "Marco", "Cd5l", "Stab2"),
  Alveolar           = c("Pparg", "Siglecf", "Itgax"),
  Microglia          = c("P2ry12", "Tmem119", "Sall1", "Hexb", "Fcrls"),
  Peritoneal_LPM     = c("Gata6", "Icam2", "F13a1"),
  Splenic_RPM        = c("Spic", "Hmox1", "Slc40a1", "Vcam1"),
  Cardiac_LYVE1_pos  = c("Lyve1", "Folr2")
)

genes_mac_curated <- lapply(genes_mac_curated, function(v) intersect(v, rownames(sce)))
genes_mac_curated <- genes_mac_curated[lengths(genes_mac_curated) > 0]

p_all_markers <- save_dotplot(sce, genes_mac_curated, "check_all_cell_markers_for_MonoMacro.pdf")
p_all_markers

assert_file("MonoMacro.sce.sub.harmony.Rds")
sce.sub <- readRDS("MonoMacro.sce.sub.harmony.Rds")

subcelltype <- data.frame(ClusterID = 0:9, subcelltype = "na")
subcelltype[subcelltype$ClusterID %in% c(0), 2]  <- "MDSC"
subcelltype[subcelltype$ClusterID %in% c(1), 2]  <- "Ly6chi Ccr2+M"
subcelltype[subcelltype$ClusterID %in% c(2), 2]  <- "APC-M"
subcelltype[subcelltype$ClusterID %in% c(3), 2]  <- "IFN response-M"
subcelltype[subcelltype$ClusterID %in% c(4), 2]  <- "Unknown"
subcelltype[subcelltype$ClusterID %in% c(5), 2]  <- "moDC"
subcelltype[subcelltype$ClusterID %in% c(6), 2]  <- "M2"
subcelltype[subcelltype$ClusterID %in% c(7), 2]  <- "inflammatory-M"
subcelltype[subcelltype$ClusterID %in% c(8), 2]  <- "Osteoclast like M"
subcelltype[subcelltype$ClusterID %in% c(9), 2]  <- "Unknown"

sce.sub <- apply_cluster_labels(sce.sub, subcelltype, "subcelltype")
write_proportion_tables(sce.sub, "subcelltype")
Idents(sce.sub) <- sce.sub$subcelltype

maincelltype <- data.frame(ClusterID = 0:9, maincelltype = "na")
maincelltype[maincelltype$ClusterID %in% c(0), 2] <- "MonoMacro"
maincelltype[maincelltype$ClusterID %in% c(1), 2] <- "MonoMacro"
maincelltype[maincelltype$ClusterID %in% c(2), 2] <- "MonoMacro"
maincelltype[maincelltype$ClusterID %in% c(3), 2] <- "MonoMacro"
maincelltype[maincelltype$ClusterID %in% c(4), 2] <- "Unknown"
maincelltype[maincelltype$ClusterID %in% c(5), 2] <- "MonoMacro"
maincelltype[maincelltype$ClusterID %in% c(6), 2] <- "MonoMacro"
maincelltype[maincelltype$ClusterID %in% c(7), 2] <- "MonoMacro"
maincelltype[maincelltype$ClusterID %in% c(8), 2] <- "MonoMacro"
maincelltype[maincelltype$ClusterID %in% c(9), 2] <- "Unknown"

sce.sub <- apply_cluster_labels(sce.sub, maincelltype, "maincelltype")
write_proportion_tables(sce.sub, "maincelltype")
Idents(sce.sub) <- sce.sub$maincelltype

sce.monomacro <- sce.sub

table(sce.monomacro$group, sce.monomacro$maincelltype)

saveRDS(sce.monomacro, "sce.MonoMacro_w_correction.Rds")
