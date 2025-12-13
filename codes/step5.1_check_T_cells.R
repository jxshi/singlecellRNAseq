# Step 5.1: check T cell subclusters

setwd("/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/3-celltype/Tcells")

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

# Load custom functions to plot.
source("~/software/functions/custom_seurat_functions.R")
source("~/software/functions/PropPlot.R")
source("~/software/functions/SubClusterPropPlot.R")
# Convert human gene symbols to mouse gene symbols.
source("~/software/functions/convertHumanGeneList.R")

# If error pops up, uncomment corresponding options here.
options(future.globals.maxSize= 891289600)
options(future.seed=TRUE)

if (!file.exists("Tcells.sce.sub.harmony.Rds")) {
  stop("Tcells.sce.sub.harmony.Rds not found in working directory")
}

sce <- readRDS("Tcells.sce.sub.harmony.Rds")
table(sce$group, sce$seurat_clusters)

# Relevel the groups in the metadata
Idents(sce) <- sce$group
sce$group <- factor(sce$group, levels = c("PBS", "DD_mGE", "DD_mIL12", "DD_mGE12", "TD_mGE12"))

# Verify the new levels
table(sce$group, sce$seurat_clusters)
Idents(sce) <- sce$seurat_clusters

dotplot_theme <- theme(
  axis.text.x = element_text(angle = 45, hjust = 1.2, vjust = 1.1, size = 8),
  strip.text  = element_text(size = 8)
)

save_dotplot <- function(object, genes, filename, width = 12, height = 10) {
  plot <- DotPlot(object, features = genes, assay = "RNA") + dotplot_theme
  ggsave(filename, plot = plot, width = width, height = height)
  plot
}

#### For mouse T cell subcluster
genes_to_check <- list(
  `γδΤ` = c("Cd3d", "Trgc2", "Trgc1"),
  `CD4_Tn` = c("Ccr7", "Sell", "Tcf7", "Il7r", "Gpr183", "Lef1", "Ltb", "Ifngr2", "S1pr1"),
  `CD4_Tcm` = c("Fas", "Hspa6", "Mt1e", "Mt1f", "Cldnd1", "Cd69", "Klf2", "Tob1", "Cxcr6", "Il6ra", "Bcl6", "Bach2"),
  `CD4_Tem` = c("Rps19", "Dusp2", "Cd44", "Gzmk"),
  `CD4_Th1` = c("Ucp2", "Arpc1b", "Tbx21", "Ifng", "Fasl", "Stmn1", "Pclaf", "Mki67", "Tnfrsf9", "Ccl5", "Cxcl13", "Pdcd1", "Icos"),
  `CD4_Th2` = c("Gata3", "Il4"),
  `CD4_Th17` = c("Linc00513", "Adam19", "Xist", "Rora", "Il17a", "Il17f", "Rorc", "Il23r", "Ctsh", "Capg", "Nr4a3"),
  `CD4_Tfh` = c( "Ctla4", "Batf", "Tnfrsf4", "Tox", "Cxcr5", "Slamf6"),
  `CD4_Treg` = c("Foxp3", "Il2ra", "Tigit", "Tnfrsf18", "Mageh1", "Sat1", "Ccr8", "Ikzf2", "Il10", "Nr4a1"),
  `CD4_Tex` = c("Prdm1", "Havcr2" ),
  `CD8_Tn` = c("Eomes", "Id3",  "S100a8", "Cst3", "Ac020916.1", "Cd28", "Cxcl10",  "Foxo1"),
  `CD8_Tcm` = c("Rps26", "Dkk3", "Btla", "Zfp36l2", "Cxcr4", "Myadm"),
  `CD8_Tem` = c( "Cx3cr1"),
  `CD8_Teff` = c("Fgfbp2", "Fcgr3a", "Klrg1", "Gzmb", "Gzma", "Gzmh", "Gzmc","Nkg7", "Satb1", "Nfkb1", "Hif1a", "Runx1", "Gnly", "Ccr2", "Ccr5", "Cxcr3", 
                 "Tnfrsf1a", "Tnfrsf1b", "Tnfsf10", "Stat1", "Isg15", "Ifit3", "Ccl3", "Ccl4", "Icos3", "Il18rap", "Il18r1", "Batf3", "Irf4",
                 "Myc", "Slc7a5", "Xcl1", "Il2"),
  `CD8_Trm组织驻留记忆` = c("Xcl2", "Tgfbr2", "Itga1", "Sipr1"),
  `CD8_Tex` = c("Lag3", "Entpd1", "Runx3", "Tox2", "Layn", "Tnfsf4", "Tim3", "Gzm", "Ung", "Mcm2", "Ccnb2", "Top2a", "Csf1", "Ccl2",
                "Tgfbi", "Cd38", "Il2rb", "Cd244", "Cd101", "Krt86", "Cd160", "Il10r", "Eome", "Nr4a2", "Ptger4", "Ifi16", "Ikzf3", "Znf683"),
  `CD8_Temra_CD45RA效应记忆` = c( "Adgrg1"),
  `CD8_细胞毒性T` = c("Prf1"),
  `CD8_TSTR应激反应` = c("Hspa1a", "Hspa1b"),
  `CD8_TSEN衰老` = c("Cd27"),
  `CD8_增值T` = c("Cdk1", "Cenpa", "Mcm5", "Pcna", "Mcm6", "Mcm3", "Tk1", "C1qb", "Mt1g", "Mt1x"),
  `CD8_活化T` = c("Cd40lg", "Anxa1"),
  `早期记忆T` = c("Zeb2")
)

p_all_markers <- DotPlot(sce, features = genes_to_check, assay = "RNA") + dotplot_theme

p_all_markers
ggsave("check_all_cell_markers_for_T.pdf", plot = p_all_markers, width = 28, height = 10)


# 3. T cells: IL7R, IL32, CD3D, CD3G, CD2, CD3E
genes_to_check <- c("IL7R", "IL32", "CD3D","CD3G","CD3E","CD2","CD4","CD8A")
genes_to_check <- convertHumanGeneList(genes_to_check)

p <- save_dotplot(sce, genes_to_check, "Tcells_featurePlot.pdf")
p



genes_to_check <- c('Sell', 'Il7r', 'Ccr7', 'Klf2', 'S1pr1', 'Tcf7',
                    'Gzmb', 'Prf1', 'Nkg7', 'Ifng', 'Il4', 'Il17a', 'Foxp3', 'Bcl6', 'Cxcr5')

p <- save_dotplot(sce, genes_to_check, "naive_central_memory_like_Tcells_featurePlot.pdf")
p




if (!require("grid")) {
  install.packages("grid")
  library(grid)
}

tmp <- table(sce$group, sce$seurat_clusters)
write.table(tmp, file = "Tcell_subtypes.csv", quote = F, sep = "\t", col.names  = NA)


# chemokines.
genes_to_check <- c("CXCL13", "CXCL12", "CCL19","CCL121","CXCL9","CXCL10","CXCL11","CCL20", "CXCL16")
genes_to_check <- convertHumanGeneList(genes_to_check)

p <- DotPlot(sce, features = genes_to_check, assay = "RNA", group.by = "group") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1.2, vjust = 1.1, size = 8),
    strip.text  = element_text(size = 8)
  )
p
ggsave(plot=p, filename="chemokines_in_Tcells.pdf",device="pdf",width = 12,height = 10)


sce.all <- readRDS("../sce.all_w_cell_type_anno.Rds")

# CD4 and CD8
genes_to_check <- c("CD4", "CD8A",'CD3D', 'CD3E', 'TRAC', 'IL7R')
genes_to_check <- convertHumanGeneList(genes_to_check)

# saveRDS(sce, "sce.Tcells.Rds")
p <- DotPlot(sce, features = genes_to_check, assay = "RNA", group.by  = "group", cols = c("lightgrey", "blue","red","green","cyan"),
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1.2, vjust = 1.1, size = 8),
    strip.text  = element_text(size = 8)
  )
p

# chemokines.
genes_to_check <- c("CXCL13", "CXCL12", "CCL19","CCL121","CXCL9","CXCL10","CXCL11","CCL20", "CXCL16")
genes_to_check <- convertHumanGeneList(genes_to_check)

p <- DotPlot(sce.all, features = genes_to_check, assay = "RNA", split.by = "group", cols = c("lightgrey", "blue","red","green","cyan"),
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1.2, vjust = 1.1, size = 8),
    strip.text  = element_text(size = 8)
  )
p
ggsave(plot=p, filename="chemokines_in_Tcells.pdf",device="pdf",width = 12,height = 10)


# chemokines.
genes_to_check <- c("Tcf7", "Tox")

p <- DotPlot(sce.all, features = genes_to_check, assay = "RNA", split.by = "group", cols = c("lightgrey", "blue","red","green","cyan"),
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1.2, vjust = 1.1, size = 8),
    strip.text  = element_text(size = 8)
  )
p
ggsave(plot=p, filename="chemokines_in_Tcells.pdf",device="pdf",width = 12,height = 10)


# chemokines.
genes_to_check <- c('Tcf1',"Tcf7", "Tox")

p <- DotPlot(sce, features = genes_to_check, assay = "RNA", split.by = "group", cols = c("lightgrey", "blue","red","green","cyan"),
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1.2, vjust = 1.1, size = 8),
    strip.text  = element_text(size = 8)
  )
p

p <- DotPlot(sce, features = genes_to_check, assay = "RNA", cols = c("lightgrey", "blue","red","green","cyan"),
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1.2, vjust = 1.1, size = 8),
    strip.text  = element_text(size = 8)
  )
p

p <- DotPlot(sce, features = genes_to_check, assay = "RNA", group.by = "group", cols = c("lightgrey", "blue","red","green","cyan"),
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1.2, vjust = 1.1, size = 8),
    strip.text  = element_text(size = 8)
  )
p

# dir.create("cyclingT")
# setwd("cyclingT")
sce.sub <- subset(sce, ident= "4")
saveRDS(sce.sub, "sce.cyclingT.Rds")

setwd("..")
dir.create("naiveTcm")
setwd("naiveTcm")
sce.sub <- subset(sce, ident= "3")
saveRDS(sce.sub, "sce.naiveTcm.Rds")

setwd("/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/3-celltype/Tcells/")
sce.sub <- readRDS("Tcells.sce.sub.harmony.Rds")
if(T) {
  table(sce.sub@active.ident)
  
  subcelltype=data.frame(ClusterID=0:15, subcelltype='na')
  
  subcelltype[subcelltype$ClusterID %in% c( 0),2]='CD4 Tn/Tcm'
  subcelltype[subcelltype$ClusterID %in% c( 1),2]='CD8 Teff'
  subcelltype[subcelltype$ClusterID %in% c( 2),2]='CD4 Teff'
  subcelltype[subcelltype$ClusterID %in% c( 3),2]='CD4 Teff'
  subcelltype[subcelltype$ClusterID %in% c( 4),2]='CD8 Teff'
  subcelltype[subcelltype$ClusterID %in% c( 5),2]='CD4 Treg'
  subcelltype[subcelltype$ClusterID %in% c( 6),2]='ProliferatingT'
  subcelltype[subcelltype$ClusterID %in% c( 7),2]='CD8 Teff'
  subcelltype[subcelltype$ClusterID %in% c( 8),2]='CD8 Teff'
  subcelltype[subcelltype$ClusterID %in% c( 9),2]='CD4 T'
  subcelltype[subcelltype$ClusterID %in% c(10),2]='CD8 Tex'
  subcelltype[subcelltype$ClusterID %in% c(11),2]='CD4 Th2'
  subcelltype[subcelltype$ClusterID %in% c(12),2]='CD4 Th17'
  subcelltype[subcelltype$ClusterID %in% c(13),2]='Unknown'
  subcelltype[subcelltype$ClusterID %in% c(14),2]='Unknown'
  subcelltype[subcelltype$ClusterID %in% c(15),2]='Unknown'
  
  
  head(subcelltype)
  subcelltype 
  table(subcelltype$subcelltype)
  sce.sub@meta.data$subcelltype = "NA"
  for(i in 1:nrow(subcelltype)){
    sce.sub@meta.data[which(sce.sub@meta.data$seurat_clusters == subcelltype$ClusterID[i]),'subcelltype'] <- subcelltype$subcelltype[i]}
  table(sce.sub@meta.data$subcelltype)
  
  cltyPerIdent <- table(sce.sub@meta.data$subcelltype, sce.sub@meta.data$orig.ident)
  write.table(cltyPerIdent, "subcelltype_per_orig.ident.tsv", quote = F, sep = "\t", row.names = T, col.names = NA)
  
  cltyPerGroup <- table(sce.sub@meta.data$subcelltype, sce.sub@meta.data$group)
  write.table(cltyPerGroup, "subcelltype_per_group.tsv", quote = F, sep = "\t", row.names = T, col.names = NA)
  
  Idents(sce.sub) <- sce.sub$subcelltype
}

if(T) {
  table(sce.sub@active.ident)
  
  maincelltype=data.frame(ClusterID=0:15, maincelltype='na')
  
  maincelltype[maincelltype$ClusterID %in% c( 0),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c( 1),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c( 2),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c( 3),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c( 4),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c( 5),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c( 6),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c( 7),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c( 8),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c( 9),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c(10),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c(11),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c(12),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c(13),2]='Unknown'
  maincelltype[maincelltype$ClusterID %in% c(14),2]='Unknown'
  maincelltype[maincelltype$ClusterID %in% c(15),2]='Unknown'
  
  
  head(maincelltype)
  maincelltype 
  table(maincelltype$maincelltype)
  sce.sub@meta.data$maincelltype = "NA"
  for(i in 1:nrow(maincelltype)){
    sce.sub@meta.data[which(sce.sub@meta.data$seurat_clusters == maincelltype$ClusterID[i]),'maincelltype'] <- maincelltype$maincelltype[i]}
  table(sce.sub@meta.data$maincelltype)
  
  cltyPerIdent <- table(sce.sub@meta.data$maincelltype, sce.sub@meta.data$orig.ident)
  write.table(cltyPerIdent, "maincelltype_per_orig.ident.tsv", quote = F, sep = "\t", row.names = T, col.names = NA)
  
  cltyPerGroup <- table(sce.sub@meta.data$maincelltype, sce.sub@meta.data$group)
  write.table(cltyPerGroup, "maincelltype_per_group.tsv", quote = F, sep = "\t", row.names = T, col.names = NA)
  
  Idents(sce.sub) <- sce.sub$maincelltype
}
sce.t <- sce.sub

table(sce.t$group,sce.t$maincelltype)

saveRDS(sce.t, "sce.T_w_correction.Rds")
