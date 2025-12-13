############### Step 8.1 Assign T cell newsubtype.

setwd("/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype/newTcells/")

setwd("/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype")
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
library(harmony)

source("~/software/functions/SubClusterPropPlot.R")
source("~/software/functions/DoHeatmapPlot.R")
source("~/software/functions/convertHumanGeneList.R")


# cname <- readRDS('celltype_name.Rds')
# fname <- readRDS('celltype_file_folder.Rds')

# for (i in 6:length(fname)) {
for (i in c(1)) {  
  setwd(fname[i])
  sce <- readRDS(paste0(fname[i], '.Rds'))
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
      FindClusters(resolution = 1.5) %>% 
      identity()
    
    umapP20 <- DimPlot(sce.sub, reduction = "umap", split.by = 'group', label = T) + plot_annotation(title = "UMAP plot of samples before harmony integration")
    tsneP20 <- DimPlot(sce.sub, reduction = "tsne", split.by = 'group', label = T) + plot_annotation(title = "TSNE plot of samples before harmony integration")
    
    library(patchwork)
    p20 <- umapP20 + tsneP20
    ggsave(filename = "umap_tsne_plots_before_using_harmony.pdf", plot = p20, device = "pdf", width = 24, height = 12)
    
    sce.sub <- sce.sub %>% RunHarmony(group.by.vars=c("group"), assay.use="RNA", reduction="pca")
    harmony_embeddings <- Embeddings(sce.sub, 'harmony')
    harmony_embeddings[1:5, 1:5]
    
    sce.sub <- sce.sub %>% 
      RunUMAP(reduction = "harmony", dims = 1:30, verbose = F) %>% 
      RunTSNE(reduction = "harmony", dims = 1:30, verbose = F) %>% 
      FindNeighbors(reduction = "harmony", k.param = 20, dims = 1:30) %>% 
      FindClusters(resolution = 1.5) %>% 
      identity()
    
    # Plot UMAP and TSNE plot after integration using harmony
    umapP2 <- DimPlot(sce.sub, reduction = "umap", split.by = 'group', label = T) + plot_annotation(title = "UMAP plot of samples after harmony integration")
    tsneP2 <- DimPlot(sce.sub, reduction = "tsne", split.by = 'group', label = T) + plot_annotation(title = "TSNE plot of samples after harmony integration")
    
    library(patchwork)
    p21 <- umapP2 + tsneP2
    ggsave(filename = "umap_tsne_plots_after_using_harmony.pdf", plot = p21, device = "pdf", width = 24, height = 12)
    
    saveRDS(sce.sub,file=paste0(fname[i],".sce.sub.harmony.Rds"))
  }
  
  # It is recommended to do differential expression on the RNA assay, and not the SCTransform.
  if(T) {
    # DefaultAssay(sce) <- "RNA"
    
    sce.markers <- FindAllMarkers(sce.sub, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
    
    sce.markers %>%
      group_by(cluster) %>%
      slice_max(n = 2, order_by = avg_log2FC)
    
    pro='_rna'
    write.csv(sce.markers,file = paste0(fname[i], pro,"_subcluster.markers.csv"))
    
    saveRDS(sce.markers, file = paste0(fname[i], pro,"_subcluster.markers.Rds"))
    
    library(dplyr)
    
    top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
    DoHeatmap(sce.sub,top10$gene,size=3)
    ggsave(filename=paste0(fname[i], pro,'_sce.markers_heatmap.pdf'),width = 14,height = 20)
    genes_to_check <- unique(top10$gene)
    # DoHeatmapPlot(sce,"batch",top10)
    p <- DotPlot(sce.sub, features = unique(top10$gene), assay='RNA')  + coord_flip()
    p
    ggsave(plot=p, filename=paste0(fname[i], pro,'DotPlot_check_top10_markers_by_clusters.pdf'), device=cairo_pdf,width = 18,height = 24)
    
    library(dplyr)
    top6 <- sce.markers %>% group_by(cluster) %>% top_n(6, avg_log2FC)
    DoHeatmap(sce.sub,top6$gene,size=3)
    ggsave(paste0(fname[i], pro,'DoHeatmap_check_top6_markers_by_clusters.pdf'),width = 18,height = 12)
    # DoHeatmapPlot(sce,"batch",top3)
    p <- DotPlot(sce.sub, features = unique(top6$gene), assay='RNA')  + coord_flip()
    p
    ggsave(plot=p, filename=paste0(fname[i], pro,'DotPlot_check_top6_markers_by_clusters.pdf'),device=cairo_pdf,width = 18,height = 12)
    
    library(dplyr)
    top3 <- sce.markers %>% group_by(cluster) %>% top_n(3, avg_log2FC)
    DoHeatmap(sce.sub,top3$gene,size=3)
    ggsave(paste0(fname[i], pro,'DoHeatmap_check_top3_markers_by_clusters.pdf'),width = 18,height = 12)
    # DoHeatmapPlot(sce,"batch",top3)
    p <- DotPlot(sce.sub, features = unique(top3$gene), assay='RNA')  + coord_flip()
    p
    ggsave(plot=p, filename=paste0(fname[i], pro,'DotPlot_check_top3_markers_by_clusters.pdf'),device=cairo_pdf,width = 18,height = 12)
    
  }
  setwd("..")
}

# Plot propotion plot on subclusters.
for (i in c(1)) {
  
  setwd(fname[i])
  sce <- readRDS(paste0(fname[i], ".sce.sub.harmony.Rds"))
  # write cell type per grouop
  cltyPerSubroup <- table(sce@meta.data$seurat_cluster, sce@meta.data$orig.ident)
  write.table(cltyPerSubroup, paste0(fname[i], "celltype_per_group.tsv"), quote = F, sep = "\t", row.names = T)
  
  # Plot proportions of different celltypes per cluster
  p <- SubPropPlot(sce, "batch") + coord_flip()
  p
  ggsave(p, filename = paste0(fname[i],"_subgroup_celltype_proportion_plot.pdf"), device = cairo_pdf, width = 10, height = 3)
  setwd("..")
}


genes_to_check <- c("Cd3d",'Cd3e',"Cd4","Cd8a",'Cd44',"Il7r",'Sell', "Lef1","Ccr7","Cxcr5","Tcf7","Gzmk","Il17a","Gata3",'Il4' )
DotPlot(sce.sub, features = genes_to_check)

# gdT and NK cell markers
Idents(sce.sub) <- sce.sub$seurat_clusters
genes_to_check <- c('Trgc1','Trgc2','Trdc','Trgv1','Trgv2',"Ncr1","Klrb1c","Itgam","Klrc1")
DotPlot(sce.sub, features = genes_to_check)

# exhausted T cell markers
genes_to_check <- c(
  "Pdcd1", "Lag3", "Havcr2", "Ctla4", "Tigit", "Cd244", "Cd160", "Entpd1",
  "Tox", "Tox2", "Nr4a1", "Nr4a2", "Eomes", "Tcf7",
  "Cxcr6", "Cxcl13", "Tnfrsf9", "Tnfrsf18", "Icos",'Cd27','Cd28','Cxcr5','Il7r'
)
DotPlot(sce.sub, features = genes_to_check)

pre_exhausted_markers <- c(
  "Tcf7",   # TCF1
  "Slamf6",
  "Cxcr5",
  "Il7r",
  "Pdcd1",
  "Tox",
  "Cd28",
  "Bcl6",
  "Btg1",
  "Btg2",
  "Ccr7",
  "Sell",
  "Lef1",
  "Id3",
  "Eomes",
  "Tnfrsf4",  # OX40
  "Tnfrsf9",  # 4-1BB
  "Icos"
)

DotPlot(sce.sub, features = pre_exhausted_markers)

################# Central Memory T cell marker genes ############
Idents(sce.sub) <- sce.sub$seurat_clusters
genes_to_check <- c("IL7R",'CCR7',"SELL","TCF7","LEF1","BCL2","GZMK","CXCR5","IL17")

genes_to_check <- convertHumanGeneList(genes_to_check)


DotPlot(sce.sub, features = genes_to_check)



sce.sub <- readRDS("Tcells.sce.sub.harmony.Rds")
if(T) {
  table(sce.sub@active.ident)
  
  newsubcelltype=data.frame(ClusterID=0:26, newsubcelltype='na')
  
  newsubcelltype[newsubcelltype$ClusterID %in% c( 0),2]='Tn'
  newsubcelltype[newsubcelltype$ClusterID %in% c( 1),2]='Cd8 Teff'
  newsubcelltype[newsubcelltype$ClusterID %in% c( 2),2]='Cd8 Tpex' #  轩??
  newsubcelltype[newsubcelltype$ClusterID %in% c( 3),2]='Th1' # Ifng, Cxcr3, Il12rb2, Tnf, Icos, Cd28, Cd40lg, Cd69, Ctla4, Entpd1
  newsubcelltype[newsubcelltype$ClusterID %in% c( 4),2]='Treg'
  newsubcelltype[newsubcelltype$ClusterID %in% c( 5),2]='Cd8 Teff'
  newsubcelltype[newsubcelltype$ClusterID %in% c( 6),2]='Cd8 Tm'
  newsubcelltype[newsubcelltype$ClusterID %in% c( 7),2]='Tfh'             # OK
  newsubcelltype[newsubcelltype$ClusterID %in% c( 8),2]='NK'                       # OK
  newsubcelltype[newsubcelltype$ClusterID %in% c( 9),2]='Tn'
  newsubcelltype[newsubcelltype$ClusterID %in% c(10),2]='Cd4 T early activated'
  newsubcelltype[newsubcelltype$ClusterID %in% c(11),2]='Proliferating Cd8 Teff'
  newsubcelltype[newsubcelltype$ClusterID %in% c(12),2]='Proliferating Cd8 T'
  newsubcelltype[newsubcelltype$ClusterID %in% c(13),2]='Cd4 Tem'
  newsubcelltype[newsubcelltype$ClusterID %in% c(14),2]='gdT'
  newsubcelltype[newsubcelltype$ClusterID %in% c(15),2]='Cd8 Tex'
  newsubcelltype[newsubcelltype$ClusterID %in% c(16),2]='Cd8 Teff'
  newsubcelltype[newsubcelltype$ClusterID %in% c(17),2]='Th2'
  newsubcelltype[newsubcelltype$ClusterID %in% c(18),2]='Proliferating Cd8 T'
  newsubcelltype[newsubcelltype$ClusterID %in% c(19),2]='Treg'
  newsubcelltype[newsubcelltype$ClusterID %in% c(20),2]='IFN-responsive Cd8 Teff'
  newsubcelltype[newsubcelltype$ClusterID %in% c(21),2]='Activated Th2'
  newsubcelltype[newsubcelltype$ClusterID %in% c(22),2]='Proliferating Treg'        # OK
  newsubcelltype[newsubcelltype$ClusterID %in% c(23),2]='Ly6c+ Cd8 Teff'
  newsubcelltype[newsubcelltype$ClusterID %in% c(24),2]='GC B'                      # OK
  newsubcelltype[newsubcelltype$ClusterID %in% c(25),2]='GC B'
  newsubcelltype[newsubcelltype$ClusterID %in% c(26),2]='Cd4 Tcm'
  

  head(newsubcelltype)
  newsubcelltype 
  table(newsubcelltype$newsubcelltype)
  sce.sub@meta.data$newsubcelltype = "NA"
  for(i in 1:nrow(newsubcelltype)){
    sce.sub@meta.data[which(sce.sub@meta.data$seurat_clusters == newsubcelltype$ClusterID[i]),'newsubcelltype'] <- newsubcelltype$newsubcelltype[i]}
  table(sce.sub@meta.data$newsubcelltype)
  
  cltyPerIdent <- table(sce.sub@meta.data$newsubcelltype, sce.sub@meta.data$orig.ident)
  write.table(cltyPerIdent, "newsubcelltype_per_orig.ident.tsv", quote = F, sep = "\t", row.names = T, col.names = NA)
  
  cltyPerGroup <- table(sce.sub@meta.data$newsubcelltype, sce.sub@meta.data$group)
  write.table(cltyPerGroup, "newsubcelltype_per_group.tsv", quote = F, sep = "\t", row.names = T, col.names = NA)
  
  Idents(sce.sub) <- sce.sub$newsubcelltype
}

if(T) {
  table(sce.sub@active.ident)
  
  maincelltype=data.frame(ClusterID=0:26, maincelltype='na')
  
  maincelltype[maincelltype$ClusterID %in% c( 0),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c( 1),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c( 2),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c( 3),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c( 4),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c( 5),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c( 6),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c( 7),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c( 8),2]='NK cells'
  maincelltype[maincelltype$ClusterID %in% c( 9),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c(10),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c(11),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c(12),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c(13),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c(14),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c(15),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c(16),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c(17),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c(18),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c(19),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c(20),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c(21),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c(22),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c(23),2]='T cells'
  maincelltype[maincelltype$ClusterID %in% c(24),2]='Unknown'
  maincelltype[maincelltype$ClusterID %in% c(25),2]='Unknown'
  maincelltype[maincelltype$ClusterID %in% c(26),2]='T cells'

  
  
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

DimPlot(sce.t, group.by = "newsubcelltype", label = T, cols = combined_colors)
DimPlot(sce.t, group.by = "newsubcelltype", label = T)
table(sce.t$newsubcelltype, sce.t$subcelltype)


saveRDS(sce.t, "sce.T_w_correction.Rds")


#################################### Plot Figures for T cells ###############################################
library(ggsci)
source("~/software/functions/groupbarcharts.R")
source("~/software/functions/groupbarchartssubcluster.R")
source("~/software/functions/PropPlot.R")
source("~/software/functions/custom_seurat_functions.R")
default_colors <- pal_nejm(palette = c("default"), alpha=1)(8)
extra_colors <- pal_aaas()(8)
combined_colors <- c(default_colors, extra_colors)
combined_colors

combined_colors <- c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#008B45FF", "#5F559BFF", "#6F99ADFF","#A20056FF",  "#FF5733",
                     "#EE4C97FF", "#3B4992FF", "#008280FF","#631879FF", "#C70039","#FFDC91FF")

combined_colors <- c("#BC3C29", "#0072B5", "#E18727", "#008B45", "#5F559B", "#6F99AD", "#A20056", "#EE4C97", "#3B4992", "#008280", "#631879", "#C70039", "#FFDC91")

combined_colors <- c(
  "#3C8A32",   # medium green
  "#28357B",  # indigo
  "#0F7374",  # teal
  "#116D69",  # dark teal
  "#5EACA7",  # turquoise
  "#B8D8DF",  # pale cyan
  "#D7E9EC",  # very light blue
  "#9BC9D9",  # blue-grey
  "#BFC7E5",  # lavender blue
  "#A9B4D7",  # grey-blue
  "#D8A5CE",  # soft pink-purple
  "#BE84B8",  # purple-pink
  "#87387D",  # plum purple
  "#D2A21B",  # mustard yellow
  "#F5C34D",  # gold
  "#C0261D",  # red
  "#A62071",  # magenta
  "#F1A800",  # bright yellow
  "#DDA873",  # peach-orange
  "#1A1B4B",  # deep navy
  "#FF5733"

)

group_colors <- c(`PBS` = "#A9B4D7",
                  `DD_mGE` = "#9BC9D9",
                  `DD_mIL12` = "#F1A800",
                  `DD_mGE12` = "#5EACA7", 
                  `TD_mGE12` = "#BE84B8")

maincelltype_colors <- c(
  `B cells`          = "#A9B4D7",  # grey-blue
  `DCs`      = "#F1A800",  # saturated golden yellow
  `NK cells`     = "#D8A5CE",  # soft pink-purple
  `T cells`       = "#BE84B8",  # purple-pink
  `Macrophages` = "#5EACA7",  # turquoise
  `Granulocytes`      = "#3C8A32",  # medium green
  `Mast cells` = "#9BC9D9"  # blue-grey
  
)

setwd("/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype/graphs/")
dir.create("TcellsFigures")

setwd("TcellsFigures/")
# sce.t <- readRDS("~/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype/sce.T_w_correction.Rds")
# sce.new <- sce.t
Idents(sce.t) <- sce.t$newsubcelltype
sce.new <- subset(sce.t, idents = setdiff(levels(Idents(sce.t)), c("GC B", "NK")))

Idents(sce.new) <- sce.new$newsubcelltype
sce.new$subcelltype <- sce.new$newsubcelltype
sce.new$group <- factor(sce.new$group, levels = c("PBS", "DD_mGE", "DD_mIL12", "DD_mGE12", "TD_mGE12"))

sce.new$subcelltype <- factor(sce.new$subcelltype, levels = c("Tn","Cd4 Tcm", "Cd4 Tem", "Th1", "Th2","Tfh", "Treg",'Proliferating Treg',
                                                              'Cd8 Tcm',"Cd8 Tem","Proliferating Cd8 T",  "Proliferating Cd8 Teff","Cd8 Teff",
                                                              "IFN-responsive Cd8 Teff","Ly6c+ Cd8 Teff", "Cd8 Tex", "gdT"))


p <- DimPlot(sce.new, group.by = "subcelltype", raster = F, label = T, cols = combined_colors)
p
ggsave("main_celltype_UMAP.pdf", plot = p, device = cairo_pdf, width = 8, height = 8)

p <- DimPlot(sce.new, group.by = "group", raster = F, label = F, cols = combined_colors)
# p <- LabelClusters(plot = p, id = "group", repel = TRUE)
p
ggsave("UMAP_by_group.pdf", plot = p, device = cairo_pdf, width = 8, height = 8)


p <- groupbarcharts(sce.new, "group")
p
ggsave("compare_maincelltype_percentage_of_each_group.pdf", plot = p, device = cairo_pdf, width = 8, height = 8)
tmp <- table(sce.new$group,sce.new$subcelltype)
write.table(tmp, file = "subcelltype_per_group.tsv", quote = F, col.names = NA, sep = "\t")

p <- groupbarchartssubcluster(sce.new, "group")
p
ggsave("compare_percentage_of_celsl_in_each_group.pdf", plot = p, device = cairo_pdf, width = 14, height = 8)

p <-  plot_integrated_clusters(sce.new, cluster_col = "subcelltype", group_col = "group")
p
ggsave("compare_percentage_of_subcelltype.pdf", plot = p, device = cairo_pdf, width = 8, height = 8)

Idents(sce.new) <- sce.new$subcelltype
DefaultAssay(sce.new) <- "RNA"
if(T) {
  # Find marker genes for each group
  sce.markers <- FindAllMarkers(sce.new, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  sce.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
  
  # DT::datatable(sce.markers)
  
  pro='cca'
  write.csv(sce.markers,file = paste0(pro,"_sce.new.markers.csv"))
  
  saveRDS(sce.markers, file = paste0(pro,"sce.new.markers.Rds"))
  
  library(dplyr) 
  top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
  DoHeatmap(sce.new,top10$gene,size=3)
  ggsave(filename=paste0(pro,'_sce.markers_heatmap.pdf'),width = 18,height = 24)
  p <- DotPlot(sce.new, features = unique(top10$gene), assay='RNA')  + coord_flip()
  p
  ggsave(plot=p, filename=paste0(pro,'DotPlot_check_top10_markers_by_clusters.pdf'), device="pdf",width = 18,height = 24)
  
  top6 <- sce.markers %>% group_by(cluster) %>% top_n(6, avg_log2FC)
  DoHeatmap(sce.new,top6$gene,size=3)
  ggsave(paste0(pro,'DoHeatmap_check_top6_markers_by_clusters.pdf'),width = 18,height = 12)
  p <- DotPlot(sce.new, features = unique(top6$gene), assay='RNA')  + coord_flip()
  p
  ggsave(plot=p, filename=paste0(pro,'DotPlot_check_top6_markers_by_clusters.pdf'),device="pdf",width = 18,height = 12)
  
  top3 <- sce.markers %>% group_by(cluster) %>% top_n(3, avg_log2FC)
  DoHeatmap(sce.new,top3$gene,size=3)
  ggsave(paste0(pro,'DoHeatmap_check_top3_markers_by_clusters.pdf'),width = 18,height = 12)
  p <- DotPlot(sce.new, features = unique(top3$gene), assay='RNA')  + coord_flip()
  p
  ggsave(plot=p, filename=paste0(pro,'DotPlot_check_top3_markers_by_clusters.pdf'),device="pdf",width = 18,height = 12)
  
  # p <- FeaturePlot(sce.new, features = c('MS4A1', 'CD79A','PTPRC','PPBP', 'HBA1', 'HBA2', 'HBB'))
  # p
  # ggsave(plot=p, filename="check_background_genes.pdf",device="pdf",width = 18,height = 12)
  
}


# progenitor-exhausted / stem-like exhausted CD8⁺ T-cell marker genes
genes_to_check <- c('Tigit', 'Ctla4', 'Pdcd1', 'Lag3', 'Havcr2', 'Cxcl13', 'Entpd1', 'Tox')
p <- DotPlot(sce.new, features = genes_to_check) + coord_flip()
p

# early memory T cell markers:
genes_to_check <- c("Cd4","Cd8a","Cd8b1","Ccr7", "Sell", "Il7r", "Tcf7", "Lef1", "Klf2")
p <- DotPlot(sce.new, features = genes_to_check) + coord_flip()
p
# VlnPlot(sce.new, features = c("Ccr7", "Sell", "Il7r", "Tcf7", "Lef1", "Klf2"), pt.size = 0)

# Effector CD8 cytotoxic T-cell signature
genes_to_check <- c("Gzmb", "Gzma", "Gzmk", "Prf1","Ifng", "Tnf","Ccl3", "Ccl4", "Ccl5","Nkg7", "Ctsw", "Cst7","Klrc1", "Klrc2",
                    "Klrd1", "Klrk1","Hopx", "Tbx21", "Zeb2")
p <- DotPlot(sce.new, features = genes_to_check) + coord_flip()
p

# Activated NK cell marker genes:
genes_to_check <- nk_markers <- c(
  "Ncr1", "Klrk1", "Klrd1", "Klrc1", "Klrc2", "Klre1",
  "Tyrobp", "Fcer1g", "Fcgr3",
  "Gzmb", "Gzma", "Gzmk", "Prf1", "Ctsw", "Cst7",
  "Ifng", "Fasl", "Tnfsf10",
  "Xcl1", "Ccl3", "Ccl4", "Ccl5",
  "Tbx21", "Eomes", "Bhlhe40", "Irf8",
  "Il2rb", "Cd160", "Klrg1", "Cd69",
  "Serpinb9", "Hilpda", "Ddit4", "Hif1a",'Mki67'
)
p <- DotPlot(sce.new, features = genes_to_check) + coord_flip()
p

# NK specific marker genes:
genes_to_check <- nk_unique_markers <- c(
  "Ncr1", "Klrb1c", 
  "Klra3", "Klra7", "Klra8", "Klra9", "Klra12", "Klra17",
  "Klrd1", "Klrc1", "Klrc2", "Klre1",
  "Tyrobp", "Fcer1g", "Fcgr3",
  "Gzmc", "Gzmg", "Gzme",
  "Eomes", "Irf8", "Cd160", "Cd244"
)
p <- DotPlot(sce.new, features = genes_to_check) + coord_flip()
p

## Treg marker genes
genes_to_check <- c("Foxp3", "Il2ra", "Ctla4", "Ikzf2", "Tnfrsf18", "Tnfrsf4","Ccr7", "Ccr4", "Ccr8", "Cxcr3",
                    "Tigit", "Havcr2", "Lag3", "Pdcd1", "Entpd1", "Nt5e","Il10", "Tgfb1", "Ebi3", "Il12a", "Nrp1","Bcl2", "Bcl2l1", "Gzmb", "Layn", "Il1rl1", "Areg")
p <- DotPlot(sce.new, features = genes_to_check) + coord_flip()
p

cd8_exhausted_markers <- c(
  "Cd8a", "Cd8b1",
  "Pdcd1", "Lag3", "Havcr2", "Cd160", "Entpd1", "Tnfrsf9",
  "Gzmf", "Gzmd", "Gzmc", "Gzmb", "Prf1", "Nkg7", "Cst7",
  "Klrd1", "Klrc1", "Klrc2", "Klrk1",
  "Ifng", "Tbx21",
  "S100a4", "Bhlhe40", "Nr4a2"
)

p <- DotPlot(sce.new, features = cd8_exhausted_markers) + coord_flip()
p

myeloid_mixed_markers <- c(
  # granulocyte / inflammatory
  "S100a8", "S100a9", "Ngp", "Retnlg", "Cxcl1", "Cxcl2",
  # monocyte / macrophage
  "Lyz2", "Cd14", "C1qa", "Aif1", "Lst1", "Apoe", "Cebpb", "Alox5ap",
  "Fcer1g", "Tyrobp", "Ctss", "Ctsd", "Acp5",
  # MHC-II
  "H2-Ab1", "H2-Aa", "H2-Eb1", "Cd74",
  # mixed T / NK features present in this cluster
  "Cd3e", "Cd8a", "Cd8b1", "Gzmk", "Nkg7", "Klrc1", "Klrc2", "Ccl5", "Cxcr3"
)
p <- DotPlot(sce.new, features = myeloid_mixed_markers) + coord_flip()
p

# Th1 and Th2 marker genes.
th1_markers <- c(
  "Tbx21",   # T-bet
  "Ifng",
  "Il2",
  "Cxcr3",
  "Ccr5",
  "Stat1",
  "Stat4",
  "Irf1",
  "Irf8",
  "Ccl3",
  "Ccl4",
  "Ccl5",
  "Gzmb",
  "Nkg7"
)

th2_markers <- c(
  "Gata3",
  "Il4",
  "Il5",
  "Il13",
  "Maf",
  "Il1rl1",  # ST2
  "Il7r",
  "Ccr4",
  "Ccl1",
  "Ccl17",
  "Ccl22",
  "Areg",
  "Csf2",
  "Ramp1",
  "Ramp3"
)

DotPlot(sce.sub, features = list(Th1 = th1_markers, Th2 = th2_markers)) + coord_flip()

sce.sub <- readRDS("~/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype/Tcells/sce.T_w_correction.Rds")
table(sce.sub$seurat_clusters)
sce.new <- sce.sub
Idents(sce.new) <- sce.new$seurat_clusters
