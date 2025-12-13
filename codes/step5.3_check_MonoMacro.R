# Step 5.2: check MonoMacrophage subclusters

setwd("/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/3-celltype/MonoMacro")

# Load the required R packages here.
library(Seurat)
library(tidyverse)
library(ggplot2)
library(sctransform) # https://satijalab.org/seurat/articles/sctransform_vignette.html
library(glmGamPoi)
library(DoubletFinder)
library(patchwork)
library(clusterProfiler)
require(org.Hs.eg.db)
# library(monocle3)
# library(garnett)
library(harmony)
library(celldex)
library(RColorBrewer)
library(patchwork)
library(future)
library(parallel)

# Load custome functions to plot.
source("~/software/functions/custom_seurat_functions.R")
source("~/software/functions/PropPlot.R")
source("~/software/functions/SubClusterPropPlot.R")
# Convert human gene symbols to mouse gene symbols.
source("~/software/functions/convertHumanGeneList.R") 

# If error pops up, uncomment corresponding options here.
options(future.globals.maxSize= 891289600)
options(future.seed=TRUE)

sce <- readRDS("MonoMacro.sce.sub.harmony.Rds")


# Curated mouse macrophage marker sets (usable with Seurat::DotPlot etc.)
genes_mac_curated <- list(
  Mono_Ly6Chi        = c("Ly6c2","Ccr2","S100a8","S100a9","Plac8","Lcn2"),
  Mono_Ly6Clo        = c("Nr4a1","Cx3cr1","Itgal","Selplg","Fcgr3"),
  
  Mac_TRM_core       = c("Adgre1","Lyz2","Mafb","Apoe","Lgals3","C1qa","C1qb","C1qc","Mertk"),
  
  M1_inflammatory    = c("Nos2","Il1b","Tnf","Il12b","Ccl2","Cxcl10","Cd86","Stat1","Irf5"),
  M2a_alternative    = c("Mrc1","Retnla","Chil3","Arg1","Ccl22","Il10ra"),
  M2c_deactivated    = c("Cd163","Tgfb1","Il10"),
  
  TAM_general        = c("Cd68","Csf1r","Maf","Axl","Vegfa","Mmp9"),
  TAM_MHCII_high     = c("H2-Aa","H2-Ab1","Cd74","Ciita","Irf1"),

  ISG_macrophage     = c("Isg15","Ifit1","Ifit3","Rsad2","Oasl1","Irf7"),
  
  LAM_Trem2_plus     = c("Trem2","Cd9","Gpnmb","Lpl","Lipa","Fabp5","Cd63","Ctsb","Ctsd"),
  SPP1_plus_TAM      = c("Spp1","Fn1"),
  
  Kupffer            = c("Clec4f","Vsig4","Timd4","Marco","Cd5l","Stab2"),
  Alveolar           = c("Pparg","Siglecf","Itgax"),
  Microglia          = c("P2ry12","Tmem119","Sall1","Hexb","Fcrls"),
  
  Peritoneal_LPM     = c("Gata6","Icam2","F13a1"),

  Splenic_RPM        = c("Spic","Hmox1","Slc40a1","Vcam1"),
  Cardiac_LYVE1_pos  = c("Lyve1","Folr2")

)

# (recommended) guard against missing features before plotting
# present <- rownames(sce)
# genes_mac_curated <- lapply(genes_mac_curated, function(v) intersect(v, present))
# genes_mac_curated <- genes_mac_curated[lengths(genes_mac_curated) > 0]


p_all_markers <- DotPlot(sce, features = genes_mac_curated, assay = "RNA") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1.2, vjust = 1.1, size = 8),
    strip.text  = element_text(size = 8)
  )

p_all_markers
ggsave("check_all_cell_markers_for_MonoMacro.pdf", plot = p_all_markers, width = 14, height = 10)



setwd("/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/3-celltype/MonoMacro")
sce.sub <- readRDS("MonoMacro.sce.sub.harmony.Rds")
if(T) {
  table(sce.sub@active.ident)
  
  subcelltype=data.frame(ClusterID=0:9, subcelltype='na')
  
  subcelltype[subcelltype$ClusterID %in% c( 0),2]='MDSC'
  subcelltype[subcelltype$ClusterID %in% c( 1),2]='Ly6chi Ccr2+M'
  subcelltype[subcelltype$ClusterID %in% c( 2),2]='APC-M'
  subcelltype[subcelltype$ClusterID %in% c( 3),2]='IFN response-M'
  subcelltype[subcelltype$ClusterID %in% c( 4),2]='Unknown'
  subcelltype[subcelltype$ClusterID %in% c( 5),2]='moDC'
  subcelltype[subcelltype$ClusterID %in% c( 6),2]='M2'
  subcelltype[subcelltype$ClusterID %in% c( 7),2]='inflammatory-M'
  subcelltype[subcelltype$ClusterID %in% c( 8),2]='Osteoclast like M'
  subcelltype[subcelltype$ClusterID %in% c( 9),2]='Unknown'

  
  
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
  Idents(sce.sub) <- sce.sub$seurat_clusters
  
  maincelltype=data.frame(ClusterID=0:9, maincelltype='na')
  
  maincelltype[maincelltype$ClusterID %in% c( 0),2]='MonoMacro'
  maincelltype[maincelltype$ClusterID %in% c( 1),2]='MonoMacro'
  maincelltype[maincelltype$ClusterID %in% c( 2),2]='MonoMacro'
  maincelltype[maincelltype$ClusterID %in% c( 3),2]='MonoMacro'
  maincelltype[maincelltype$ClusterID %in% c( 4),2]='Unknown'
  maincelltype[maincelltype$ClusterID %in% c( 5),2]='MonoMacro'
  maincelltype[maincelltype$ClusterID %in% c( 6),2]='MonoMacro'
  maincelltype[maincelltype$ClusterID %in% c( 7),2]='MonoMacro'
  maincelltype[maincelltype$ClusterID %in% c( 8),2]='MonoMacro'
  maincelltype[maincelltype$ClusterID %in% c( 9),2]='Unknown'

  
  
  
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
sce.monomacro <- sce.sub

table(sce.monomacro$group,sce.monomacro$maincelltype)

saveRDS(sce.monomacro, "sce.MonoMacro_w_correction.Rds")

