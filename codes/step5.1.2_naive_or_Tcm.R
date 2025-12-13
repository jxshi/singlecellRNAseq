#################### Step 5.1.2. naive_or_Tcm cells ####################
library(Seurat)
library(tidyverse)
library(reticulate)
library(ggplot2)
library(patchwork)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)

source("~/software/functions/SubClusterPropPlot.R")
source("~/software/functions/DoHeatmapPlot.R")

setwd("/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/3-celltype/Tcells/naiveTcm")

sce <- readRDS("sce.naiveTcm.Rds")

DefaultAssay(sce) <- "RNA"

sceList <- SplitObject(sce, split.by = "group")

if(T) {
  # Batch correction using harmony
  for (j in seq(length(sceList))) {
    DefaultAssay(sceList[[j]]) <- "RNA"
    sceList[[j]] <- NormalizeData(sceList[[j]]) %>% 
      FindVariableFeatures(selection.method = "vst", nfeatures = 3000)
  }
  sce.sub <- merge(sceList[[1]], y= sceList[-1]) 
  DefaultAssay(sce.sub) <- "RNA"
  
  sce.sub <- NormalizeData(sce.sub) %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30, umap.method = "uwot") %>%
    RunTSNE(reduction = "pca", dims = 1:30) %>%
    FindNeighbors(reduction = "pca", k.param = 20, dims = 1:30) %>% 
    FindClusters(resolution = 0.8) %>% 
    identity()
  
  umapP20 <- DimPlot(sce.sub, reduction = "umap", split.by = 'group', label = T) + plot_annotation(title = "UMAP plot of samples before harmony integration")
  tsneP20 <- DimPlot(sce.sub, reduction = "tsne", split.by = 'group', label = T) + plot_annotation(title = "TSNE plot of samples before harmony integration")
  
  library(patchwork)
  p20 <- umapP20 + tsneP20
  ggsave(filename = "umap_tsne_plots_before_using_harmony.pdf", plot = p20, device = "pdf", width = 24, height = 12)
  
  sce.sub <- sce.sub %>% RunHarmony(group.by.vars=c("orig.ident","group"), assay.use="RNA", reduction="pca")
  harmony_embeddings <- Embeddings(sce.sub, 'harmony')
  harmony_embeddings[1:5, 1:5]
  
  sce.sub <- sce.sub %>% 
    RunUMAP(reduction = "harmony", dims = 1:30, verbose = F) %>% 
    RunTSNE(reduction = "harmony", dims = 1:30, verbose = F) %>% 
    FindNeighbors(reduction = "harmony", k.param = 20, dims = 1:30) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()
  
  # Plot UMAP and TSNE plot after integration using harmony
  umapP2 <- DimPlot(sce.sub, reduction = "umap", split.by = 'group', label = T) + plot_annotation(title = "UMAP plot of samples after harmony integration")
  tsneP2 <- DimPlot(sce.sub, reduction = "tsne", split.by = 'group', label = T) + plot_annotation(title = "TSNE plot of samples after harmony integration")
  
  library(patchwork)
  p21 <- umapP2 + tsneP2
  ggsave(filename = "umap_tsne_plots_after_using_harmony.pdf", plot = p21, device = "pdf", width = 24, height = 12)
  
  saveRDS(sce.sub,file=paste0("naiveTcm",".sce.sub.harmony.Rds"))
}

# It is recommended to do differential expression on the RNA assay, and not the SCTransform.
if(T) {
  DefaultAssay(sce) <- "RNA"
  
  sce.markers <- FindAllMarkers(sce.sub, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
  
  sce.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
  
  pro='_rna'
  write.csv(sce.markers,file = paste0('naiveTcm', pro,"_subcluster.markers.csv"))
  
  saveRDS(sce.markers, file = paste0('naiveTcm', pro,"_subcluster.markers.Rds"))
  
  library(dplyr)
  
  top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
  DoHeatmap(sce.sub,top10$gene,size=3)
  ggsave(filename=paste0('naiveTcm', pro,'_sce.markers_heatmap.pdf'),width = 14,height = 20)
  genes_to_check <- unique(top10$gene)
  # DoHeatmapPlot(sce,"batch",top10)
  p <- DotPlot(sce.sub, features = unique(top10$gene), assay='RNA')  + coord_flip()
  p
  ggsave(plot=p, filename=paste0('naiveTcm', pro,'DotPlot_check_top10_markers_by_clusters.pdf'), device=cairo_pdf,width = 18,height = 24)
  
  library(dplyr)
  top6 <- sce.markers %>% group_by(cluster) %>% top_n(6, avg_log2FC)
  DoHeatmap(sce.sub,top6$gene,size=3)
  ggsave(paste0('naiveTcm', pro,'DoHeatmap_check_top6_markers_by_clusters.pdf'),width = 18,height = 12)
  # DoHeatmapPlot(sce,"batch",top3)
  p <- DotPlot(sce.sub, features = unique(top6$gene), assay='RNA')  + coord_flip()
  p
  ggsave(plot=p, filename=paste0('naiveTcm', pro,'DotPlot_check_top6_markers_by_clusters.pdf'),device=cairo_pdf,width = 18,height = 12)
  
  library(dplyr)
  top3 <- sce.markers %>% group_by(cluster) %>% top_n(3, avg_log2FC)
  DoHeatmap(sce.sub,top3$gene,size=3)
  ggsave(paste0('naiveTcm', pro,'DoHeatmap_check_top3_markers_by_clusters.pdf'),width = 18,height = 12)
  # DoHeatmapPlot(sce,"batch",top3)
  p <- DotPlot(sce.sub, features = unique(top3$gene), assay='RNA')  + coord_flip()
  p
  ggsave(plot=p, filename=paste0('naiveTcm', pro,'DotPlot_check_top3_markers_by_clusters.pdf'),device=cairo_pdf,width = 18,height = 12)
  
}

# CD4 and CD8
genes_to_check <- c("CD4", "CD8A",'CD3D', 'CD3E', 'TRAC', 'IL7R','CD44','IL2RA','FOXP3',"CTLA4",'GZMB', 'PRF1', 'CCL5', 'NKG7', 'GZMA')
genes_to_check <- convertHumanGeneList(genes_to_check)

# saveRDS(sce, "sce.Tcells.Rds")
p <- DotPlot(sce.sub, features = genes_to_check, assay = "RNA",  cols = c("lightgrey", "blue","red","green","cyan"),
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1.2, vjust = 1.1, size = 8),
    strip.text  = element_text(size = 8)
  )
p
ggsave(p, filename = paste0('naiveTcm',"_Dotplot_plot.pdf"), device = cairo_pdf, width = 10, height = 10)

#Key Marker Genes for γδ T Cells
genes_to_check <- c('Tcrg-C1', 'Tcrg-C2', 'Tcrg-C3', 'Tcrg-C4', 'Trac', 'Trbc1', 'Trbc2')
p <- DotPlot(sce.sub, features = genes_to_check, assay = "RNA",  cols = c("lightgrey", "blue","red","green","cyan"),
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1.2, vjust = 1.1, size = 8),
    strip.text  = element_text(size = 8)
  )
p

if(T) {
  table(sce.sub@active.ident)
  
  subcelltype=data.frame(ClusterID=0:10, subcelltype='na')
  
  subcelltype[subcelltype$ClusterID %in% c( 0),2]='CD4_Tcm'
  subcelltype[subcelltype$ClusterID %in% c( 1),2]='CD8_Tem'
  subcelltype[subcelltype$ClusterID %in% c( 2),2]='CD8_Tem'
  subcelltype[subcelltype$ClusterID %in% c( 3),2]='CD4_Tem'
  subcelltype[subcelltype$ClusterID %in% c( 4),2]='Tem'
  subcelltype[subcelltype$ClusterID %in% c( 5),2]='eTreg'
  subcelltype[subcelltype$ClusterID %in% c( 6),2]='CD4_Tcm'
  subcelltype[subcelltype$ClusterID %in% c( 7),2]='CD4_Tcm'
  subcelltype[subcelltype$ClusterID %in% c( 8),2]='CD8_Tem'
  subcelltype[subcelltype$ClusterID %in% c( 9),2]='CD8_Teff'
  subcelltype[subcelltype$ClusterID %in% c(10),2]='CD8_Tem'

  
  
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
  
  saveRDS(sce.sub, "naiveTcm_w_subcelltype.Rds")
}



# Plot propotion plot on subclusters.
if(T) {
  
  sce <- readRDS(paste0('naiveTcm', ".sce.sub.harmony.Rds"))
  # write cell type per grouop
  cltyPerSubroup <- table(sce@meta.data$seurat_cluster, sce@meta.data$orig.ident)
  write.table(cltyPerSubroup, paste0('naiveTcm', "celltype_per_group.tsv"), quote = F, sep = "\t", row.names = T)
  
  # Plot proportions of different celltypes per cluster
  p <- SubPropPlot(sce, "batch") + coord_flip()
  p
  ggsave(p, filename = paste0('naiveTcm',"_subgroup_celltype_proportion_plot.pdf"), device = cairo_pdf, width = 10, height = 3)
}

# For MonoMacro subcluster
# GeneMarkers from Li Wei
if(T) {
  # All interested gene markers in Mono/Macro cluster
  genes_to_check = c('CD14', 'PTPRC', 'CD68', # Universal markers
                     'IRF1','IRF2BP2','IFI27L2','IFI30','IFITM3','IRF5','IFI44L','MX1','IFIT1','IFI44','ISG15','IFIT3','IFIT2', 'IFI6',
                     'IFI16','IFI35','IFIH1','IRF7', 'IRF8','IRF3','IFITM2','ISG20L2', # Interferon related
                     'HLA-A', 'HLA-B','HLA-C','HLA-E','HLA-F', # MHC-I
                     'HLA-DRB1', 'HLA-DRA', 'HLA-DQB1', 'HLA-DPB1', 'HLA-DPA1', 'CD74', 'HLA-DMA', 'HLA-DMB','HLA-DRB5', # MHC-II
                     'FCGR3A', 'FCGR1A', 'FCGR3B', 'FCGR2A','FCGR2B','FCGR2C',  # Fc receptor genes
                     'CD163', 'CD36', 'MARCO','MRC1', # Savenger receptor
                     'SIRPA', 'SIGLEC10', 'LILRB1', 'LILRB2', 'PDCD1', 'SLAMF7', # Don't eat me markers
                     'MPP7','MMP9' ,'MMP8', # MMP genes
                     'S100A4', 'S100A6', 'S100A8', 'S100A9','S100A10','S100A11', 'S100A12', # S100 genes
                     'GAS5', 'GAS7', 'VEGFA','VEGFB' , 'IL1B', 'IL18', # cytokines
                     'MERTK','AXL','TIMD4','SLC40A1', 'HMOX1', #吞噬和铁循环
                     'ANXA1','HMGB2', 'NLRP3', 'AIF1', 'SIGLEC1', 'RNASE2', #other markers
                     'TPM2','LILRA4', # Dendritic cells
                     'LAMP3', 'IDO1','IDO2',## DC3
                     'CD1E','CD1C') # DC2)
  
  p <- DotPlot(sce, features = genes_to_check, assay='RNA')  + coord_flip()
  p
  ggsave(plot=p, filename=paste0('naiveTcm', pro,'DotPlot_check_interested_markers_by_clusters_MonoMacro_cells.pdf'), device=cairo_pdf,width = 16,height = 24)
  
  # Universal and MHC-I MHC-II gene markers
  genes_to_check = c('CD14', 'PTPRC', 'CD68', # Universal markers
                     'HLA-A', 'HLA-B','HLA-C','HLA-E','HLA-F', # MHC-I
                     'HLA-DRB1', 'HLA-DRA', 'HLA-DQB1', 'HLA-DPB1', 'HLA-DPA1', 'CD74', 'HLA-DMA', 'HLA-DMB','HLA-DRB5') # MHC-II
  
  p <- DotPlot(sce, features = genes_to_check, assay='RNA')  + coord_flip() + ggtitle('Universal makers and MHC gene markers')
  p
  ggsave(plot=p, filename=paste0('naiveTcm', pro,'DotPlot_univ_mhc_by_clusters.pdf'), device=cairo_pdf,width=9,height=10)
  
  # Interferon gene markers
  genes_to_check = c('IRF1','IRF2BP2','IFI27L2','IFI30','IFITM3','IRF5','IFI44L','MX1','IFIT1','IFI44','ISG15','IFIT3','IFIT2', 'IFI6',
                     'IFI16','IFI35','IFIH1','IRF7', 'IRF8','IRF3','IFITM2','ISG20L2') # Interferon related
  
  p <- DotPlot(sce, features = genes_to_check, assay='RNA')  + coord_flip() + ggtitle('Interferon genes')
  p
  ggsave(plot=p, filename=paste0('naiveTcm', pro,'DotPlot_interferon_markers_by_clusters.pdf'), device=cairo_pdf,width=9,height=10)
  
  
  genes_to_check = c('FCGR3A', 'FCGR1A', 'FCGR3B', 'FCGR2A','FCGR2B','FCGR2C',  # Fc receptors
                     'CD163', 'CD36', 'MARCO','MRC1', # Scavenger receptors
                     'SIRPA', 'SIGLEC10', 'LILRB1', 'LILRB2', 'PDCD1', 'SLAMF7', # Don't eat me markers
                     'MPP7','MMP9' ,'MMP8') # MMP genes
  
  p <- DotPlot(sce, features = genes_to_check, assay='RNA')  + coord_flip() + ggtitle('Fc Receptors, Scavenger Receptors, Don\'t eat me and MMP genes')
  p
  ggsave(plot=p, filename=paste0('naiveTcm', pro,'DotPlot_fcReceptor_markers_by_clusters.pdf'), device=cairo_pdf,width=9,height=10)
  
  genes_to_check = c('S100A4', 'S100A6', 'S100A8', 'S100A9','S100A10','S100A11', 'S100A12', # S100 genes
                     'GAS5', 'GAS7', 'VEGFA','VEGFB' , 'IL1B', 'IL18', # cytokines
                     'MERTK','AXL','TIMD4','SLC40A1', 'HMOX1', #吞噬和铁循环
                     'ANXA1','HMGB2', 'NLRP3', 'AIF1', 'SIGLEC1', 'RNASE2') #other markers
  
  p <- DotPlot(sce, features = genes_to_check, assay='RNA')  + coord_flip() + ggtitle('S100 gene family, cytokines and other markers')
  p
  ggsave(plot=p, filename='DotPlot_s100_cytokine_markers_by_clusters.pdf', device=cairo_pdf,width=9,height=10)
}

### For B and PlasmaCells
# GeneMarkers from Li Wei
setwd("Plasmacells")
sce <- readRDS("Plasmacellssubcluster_best_resolution.Rds")
# CD20=MS4A1; CD138=SDC1; BCMA=TNFRSF17; IGL3=IGLC3?
genes_to_check = c('MS4A1','CD79B','CD79A','CD22','IGHM','IGHD','IGLC3','CD19','CD38','SDC1','TNFRSF17','SLAMF7','GPRC5D','CD24',
                   'GAS6','GAS7','IGF1','CD27')
p <- DotPlot(sce, features = genes_to_check, assay='SCT')  + coord_flip()
p

ggsave(plot=p, filename='DotPlot_check_interested_markers_for_B_Plasma_cells.pdf', device=cairo_pdf,width = 10,height = 8)
