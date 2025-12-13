#################### Step 3.1 Cell Type Recognition ####################
if(!dir.exists("../3-celltype")){
  dir.create("../3-celltype")
}
setwd("../3-celltype")

library(celldex)
library(SingleR)
library(Seurat)
library(gplots)
library(tidyverse)

sce.all <- readRDS("../1-QC/sce.all_best_resolution.Rds")

# Use singleR to annotate cell type to help identify cell types
if(T){
  
  sce <- as.SingleCellExperiment(DietSeurat(sce.all))
  sce
  
  # Choose from monaco, hpca, dice, or use all of them.
  anno <- c("monaco", "hpca", "dice")
  
  if("monaco" %in% anno) {
    monaco.ref <- celldex::MonacoImmuneData()
    
    monaco.main <- SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.main)
    monaco.fine <- SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.fine)
    
    table(monaco.main$pruned.labels)
    table(monaco.fine$pruned.labels)
    
    # Add the annotations to the Seurat object metadata
    sce.all@meta.data$monaco.main <- monaco.main$pruned.labels
    sce.all@meta.data$monaco.fine <- monaco.fine$pruned.labels
    
    monaco.main.tab <- table(sce.all@meta.data$monaco.main, sce.all@meta.data$seurat_clusters)
    monaco.fine.tab <- table(sce.all@meta.data$monaco.fine, sce.all@meta.data$seurat_clusters)
    
    write.table(monaco.main.tab, file = "monaco.main.table_for_sce.all.tsv", sep="\t", quote = F)
    write.table(monaco.fine.tab, file = "monaco.fine.table_for_sce.all.tsv", sep="\t", quote = F)
    
    pdf( "monaco.main.label_balloonplot.pdf", width = 16, height = 10)
    balloonplot(monaco.main.tab)
    dev.off()
    pdf( "monaco.fine.label_balloonplot.pdf", width = 16, height = 10)
    alloonplot(monaco.fine.tab)
    dev.off()
    
  }
  
  if("hpca" %in% anno) {
    hpca.ref <- celldex::HumanPrimaryCellAtlasData()
    
    hpca.main <- SingleR(test = sce, assay.type.test = 1, ref = hpca.ref, labels = hpca.ref$label.main)
    hpca.fine <- SingleR(test = sce, assay.type.test = 1, ref = hpca.ref, labels = hpca.ref$label.fine)
    
    table(hpca.main$pruned.labels)
    table(hpca.fine$pruned.labels)
    
    # Add the annotations to the Seurat object metadata
    sce.all@meta.data$hpca.main   <- hpca.main$pruned.labels
    sce.all@meta.data$hpca.fine   <- hpca.fine$pruned.labels
    
    hpca.main.tab <- table(sce.all@meta.data$hpca.main, sce.all@meta.data$seurat_clusters)
    hpca.fine.tab <- table(sce.all@meta.data$hpca.fine, sce.all@meta.data$seurat_clusters)
    
    write.table(hpca.main.tab, file = "hpca.main.table_for_sce.all.tsv", sep="\t", quote = F)
    write.table(hpca.fine.tab, file = "hpca.fine.table_for_sce.all.tsv", sep="\t", quote = F)
    
    pdf( "hpca.main.label_balloonplot.pdf", width = 16, height = 10)
    balloonplot(hpca.main.tab)
    dev.off()
    pdf( "hpca.fine.label_balloonplot.pdf", width = 16, height = 10)
    balloonplot(hpca.fine.tab)
    dev.off()
  }
  
  if("dice" %in% anno) {
    dice.ref <- celldex::DatabaseImmuneCellExpressionData()
    
    dice.main <- SingleR(test = sce, assay.type.test = 1, ref = dice.ref, labels = dice.ref$label.main)
    dice.fine <- SingleR(test = sce, assay.type.test = 1, ref = dice.ref, labels = dice.ref$label.fine)
    
    table(dice.main$pruned.labels)
    table(dice.fine$pruned.labels)
    
    # Add the annotations to the Seurat object metadata
    sce.all@meta.data$dice.main   <- dice.main$pruned.labels
    sce.all@meta.data$dice.fine   <- dice.fine$pruned.labels
    
    dice.main.tab <- table(sce.all@meta.data$dice.main, sce.all@meta.data$seurat_clusters)
    dice.fine.tab <- table(sce.all@meta.data$dice.fine, sce.all@meta.data$seurat_clusters)
    
    write.table(dice.main.tab, file = "dice.main.table_for_sce.all.tsv", sep="\t", quote = F)
    write.table(dice.fine.tab, file = "dice.fine.table_for_sce.all.tsv", sep="\t", quote = F)
    
    pdf( "dice.main.label_balloonplot.pdf", width = 16, height = 10)
    balloonplot(dice.main.tab)
    dev.off()
    pdf( "dice.fine.label_balloonplot.pdf", width = 16, height = 10)
    balloonplot(dice.fine.tab)
    dev.off()
  }
}

#################### Step 3.2 Assign cell type to each cluster ####################
# **************** IMPORTANT ****************
#This part has to be done manually.
# Option 1:
if(T) {
  table(sce.all@active.ident)
  
  celltype=data.frame(ClusterID=0:28, celltype='na')
  
  celltype[celltype$ClusterID %in% c( 0),2]='T cells'
  celltype[celltype$ClusterID %in% c( 1),2]='T cells'
  celltype[celltype$ClusterID %in% c( 2),2]='T cells'
  celltype[celltype$ClusterID %in% c( 3),2]='T cells'
  celltype[celltype$ClusterID %in% c( 4),2]='Mono Macro'
  celltype[celltype$ClusterID %in% c( 5),2]='Mono Macro'
  celltype[celltype$ClusterID %in% c( 6),2]='Mono Macro'
  celltype[celltype$ClusterID %in% c( 7),2]='Mono Macro'
  celltype[celltype$ClusterID %in% c( 8),2]='Mono Macro'
  celltype[celltype$ClusterID %in% c( 9),2]='B cells'
  celltype[celltype$ClusterID %in% c(10),2]='T cells'
  celltype[celltype$ClusterID %in% c(11),2]='T cells'
  celltype[celltype$ClusterID %in% c(12),2]='T cells'
  celltype[celltype$ClusterID %in% c(13),2]='DCs'
  celltype[celltype$ClusterID %in% c(14),2]='Mono Macro'
  celltype[celltype$ClusterID %in% c(15),2]='T cells'
  celltype[celltype$ClusterID %in% c(16),2]='T cells'
  celltype[celltype$ClusterID %in% c(17),2]='Mono Macro'
  celltype[celltype$ClusterID %in% c(18),2]='Mono Macro'
  celltype[celltype$ClusterID %in% c(19),2]='Epithelial cells'
  celltype[celltype$ClusterID %in% c(20),2]='T cells'
  celltype[celltype$ClusterID %in% c(21),2]='T cells'
  celltype[celltype$ClusterID %in% c(22),2]='Mast cells'
  celltype[celltype$ClusterID %in% c(23),2]='T cells'
  celltype[celltype$ClusterID %in% c(24),2]='DCs'
  celltype[celltype$ClusterID %in% c(25),2]='DCs'
  celltype[celltype$ClusterID %in% c(26),2]='Mono Macro'
  celltype[celltype$ClusterID %in% c(27),2]='Fibroblasts'
  celltype[celltype$ClusterID %in% c(28),2]='Neutrophils'
  
  head(celltype)
  celltype 
  table(celltype$celltype)
  sce.all@meta.data$celltype = "NA"
  for(i in 1:nrow(celltype)){
    sce.all@meta.data[which(sce.all@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
  table(sce.all@meta.data$celltype)
  
  cltyPerIdent <- table(sce.all@meta.data$celltype, sce.all@meta.data$orig.ident)
  write.table(cltyPerIdent, "celltype_per_orig.ident.tsv", quote = F, sep = "\t", row.names = T, col.names = NA)
  
  cltyPerGroup <- table(sce.all@meta.data$celltype, sce.all@meta.data$group)
  write.table(cltyPerGroup, "celltype_per_group.tsv", quote = F, sep = "\t", row.names = T, col.names = NA)
  
  Idents(sce.all) <- sce.all$celltype
}

# Plot proportion plot of cell types
if(T) {
  source("~/software/functions/PropPlot.R")
  # Plot proportions of different celltypes per group
  p <- PropPlot(sce.all, "group") # + coord_flip()
  p1 <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(p1, filename = "group_celltype_proportion_plot.pdf", device = cairo_pdf, width = 6, height = 10)
  
  # Plot proportions of different celltypes per orig.ident
  p <- PropPlot(sce.all, "batch")
  p1 <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(p1, filename = "batch_celltype_proportion_plot.pdf", device = cairo_pdf, width = 6, height = 10)
}

# Option 2:
if(T) {
  sce.all<- RenameIdents(sce.all,
                         `0` = "Mono/Macro",
                         `1` = "B cells", #done
                         `2` = "CD4+ Memory T", #done
                         `3` = "CD8+ NK cells", #done
                         `4` = "NK cells", # ?
                         `5` = "Fibroblasts", # MME?
                         `6` = "Endothelial cells", #done PECAM1 ITGAX
                         `7` = "Mono/Macro cells", #done
                         `8` = "Fibroblasts", # MME ?
                         `9` = "Myeloid", #MAF
                         `10` = "Mono/Macro cells", #done
                         `11` = "Myeloid", # ITGAE
                         `12` = "Naive B cells", #done TCL1A MME
                         `13` = "Plasma B cells", #done TNFRSF17
                         `14` = "DC", # CD1E, CD1C, ITGAX
                         `15` = "NK cells") #? GZMB
  DimPlot(sce.all, label = TRUE)
  saveRDS(sce.all,"sce.all_w_cell_type_anno.Rds")
}
# **************** IMPORTANT ****************
#################### Step 3.3 Check known markers for assigned types of cells ####################
if(T) {
  library(ggplot2) 
  genes_to_check = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A','CD19', 'CD79A', 'MS4A1' ,
                     'IGHG1', 'MZB1', 'SDC1',
                     'CD68', 'CD163', 'CD14','CD206',
                     'TPSAB1' , 'TPSB2',  # mast cells,
                     'RCVRN','FPR1' , 'ITGAM' ,
                     'C1QA',  'C1QB',  # mac
                     'S100A9', 'S100A8', 'MMP19',# monocyte
                     'LAMP3', 'IDO1','IDO2','ITGAX','CD74', 'FLT3','CLEC9A','CLEC10A','TCF4',## DC3 
                     'CD1E','CD1C', # DC2
                     'KLRB1','NCR1', # NK 
                     'FGF7','MME', 'ACTA2',
                     'PECAM1', 'VWF', 
                     'MKI67','TOP2A',
                     'EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1' )
  genes_to_check <- convertHumanGeneList(genes_to_check)
  library(stringr)   
  p_all_markers <- DotPlot(sce.all, features = unique(genes_to_check),
                           group.by = "seurat_clusters",assay='RNA')  + coord_flip()
  
  p_all_markers
  ggsave(plot=p_all_markers, filename="check_all_marker_by_seurat_clusters.pdf", width = 14, height = 10) 
  
  p <- DimPlot(sce.all, reduction = "umap", group.by = "celltype",label = T, raster=FALSE)
  p
  ggsave(plot = p, filename='umap_by_celltype.pdf', width = 14, height = 10)
  
  tab.1=table(sce.all@meta.data$celltype,sce.all@meta.data$seurat_clusters) 
  library(gplots)
  tab.1
  pro='cluster'
  pdf(file = paste0(pro,'_celltype_vs_seurat_clusters.pdf'), width = 14, height = 10)
  balloonplot(tab.1, main =" Celltype VS Seurat clusters ", xlab ="", ylab="",
              label = T, show.margins = F)
  dev.off()
  
  library(patchwork)
  th=theme(axis.text.x = element_text(angle = 45, 
                                      vjust = 0.5, hjust=0.5))
  p_all_markers=DotPlot(sce.all, features = unique(genes_to_check),
                        assay='RNA' ,group.by = 'celltype' )  + coord_flip()+ th
  p_umap <- DimPlot(sce.all, reduction = "umap", group.by = "celltype",label = T, raster=FALSE) 
  p <- p_all_markers+p_umap
  p
  ggsave(plot = p, filename = 'markers_umap_by_celltype.pdf', device = "pdf", width = 16, height = 12)
  
  p <- DimPlot(sce.all, reduction = "umap", group.by = "celltype", split.by = "orig.ident", label = T, raster=FALSE) 
  p
  ggsave(plot = p, filename = 'umap_by_celltype_split_orig.ident.pdf', width = 24, height = 12)
  
  saveRDS(sce.all, "sce.all_w_cell_type_anno.Rds")
}

if(T) {
  # sce.all <- readRDS("sce.all_w_cell_type_anno.Rds")
  phe=sce.all@meta.data
  saveRDS(phe,file = 'phe_by_markers.Rds')
  
  sce=sce.all
  Idents(sce)=sce$celltype
  DefaultAssay(sce) <- "RNA"
  table(Idents(sce)) 
  sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
  DT::datatable(sce.markers)
  
  pro='celltype_deg'
  
  write.csv(sce.markers,file=paste0(pro,'_sce.markers.csv'))
  
  library(dplyr) 
  # Plot heatmap and dotplot for top 10 markers of each cluster.
  top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
  p <- DoHeatmap(sce,top10$gene,size=3)
  p
  ggsave(plot = p, filename=paste0(pro,'_sce.markers_check_top10_heatmap.pdf'), device = "pdf",width = 14, height = 10)
  p <- DotPlot(sce, features = unique(top10$gene),
               assay='RNA')  + coord_flip()
  p
  ggsave(plot = p, paste0(pro,'DotPlot_check_top10_markers_by_clusters.pdf'), device = "pdf", width = 14, height = 10)
  
  library(dplyr) 
  # Plot heatmap and dotplot for top 6 markers of each cluster.
  top6 <- sce.markers %>% group_by(cluster) %>% top_n(6, avg_log2FC)
  p <- DoHeatmap(sce,top6$gene,size=3)
  ggsave(p, paste0(pro,'DoHeatmap_check_top6_markers_by_clusters.pdf'),width = 18,height = 12)
  # DoHeatmapPlot(sce,"batch",top3)
  p <- DotPlot(sce, features = unique(top6$gene), assay='RNA')  + coord_flip()
  p
  ggsave(plot=p, filename=paste0(pro,'DotPlot_check_top6_markers_by_clusters.pdf'),device=cairo_pdf,width = 18,height = 12)
  
  
  library(dplyr) 
  # Plot heatmap and dotplot for top 3 markers of each cluster.
  top3 <- sce.markers %>% group_by(cluster) %>% top_n(3, avg_log2FC)
  p <- DoHeatmap(sce,top3$gene,size=3)
  p
  ggsave(paste0(pro,'DoHeatmap_check_top3_markers_by_clusters.pdf'), device = "pdf", width = 14, height = 10)
  p <- DotPlot(sce, features = unique(top3$gene),
               assay='RNA')  + coord_flip()
  p
  ggsave(plot = p, filename = paste0(pro,'DotPlot_check_top3_markers_by_clusters.pdf'), device = "pdf", width = 14, height = 10)
  saveRDS(sce.markers,file = paste0(pro,'sce.cell.type.markers.Rds'))
}

# Feaureplots of interested marker genes
if(T) {
  DefaultAssay(sce) <- "RNA"
  genes_to_check <- c('CD3D','CD4','CD8A','GNLY','IL7R','GZMK','FCGR3A','CD14','CD68','CD163','TNFRSF17','CD19',
                      'CD22','MS4A1','CD79A','IGKC','CD1C','TPM2')
  genes_to_check <- convertHumanGeneList(genes_to_check)
  
  p <- FeaturePlot(sce, genes_to_check, raster=FALSE)
  p
  ggsave(plot = p, filename = "major_cluster_feature_plot.pdf",device = cairo_pdf, width = 14, height = 14)
}

#################### Step 3.4 Split data according to cell types ####################
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)

#拆分为 多个 seurat子对象
sce.all.list <- SplitObject(sce.all, split.by = "celltype")
sce.all.list
cname <- names(sce.all.list)
fname <- gsub(' ', '', names(sce.all.list))
fname
saveRDS(cname, "celltype_name.Rds")
saveRDS(fname, "celltype_file_folder.Rds")

for (i in fname) {
  dir.create(i)
}

for (i in seq(length(sce.all.list))) {
  sce.sub <- paste0("sce.",fname[i])
  sce.sub <- subset(sce.all, idents =cname[i])
  sce.sub
  table(sce.sub$celltype, sce.sub$group) 
  saveRDS(sce.sub,file = paste0(fname[i],'/',fname[i],'.Rds'))
}


### DotPlot for selected genes

interested_markers <- c("CD14", "CD68", "LYZ", "FCGR3B", "CXCR2", "S100A8","S100A9","IL7R", "IL32", "CD3D","CD3G","CD3E","CD2",
                        "IGHM", "MS4A1", "CD79A","CD79B","IGHD","GNLY", "NKG7", "GZMB","KLRD1","FCER1A", "CD1C", "MRC1",
                        "GATA1", "GATA2")
library(stringr)   

# Plot by original cluster ID
p_interested_markers <- DotPlot(sce.all, features = unique(interested_markers),
                                group.by = "seurat_clusters",assay='RNA')  + coord_flip()

p_interested_markers
ggsave(plot=p_interested_markers, filename="check_interested_marker_by_seurat_clusters.pdf", width = 14, height = 10) 

# plot by assigned cell type
p_interested_markers <- DotPlot(sce.all, features = unique(interested_markers),
                                group.by = "celltype",assay='RNA')  + coord_flip()

p_interested_markers
ggsave(plot=p_interested_markers, filename="check_interested_marker_by_assigned_cell_type.pdf", width = 10, height = 6) 

p <- DimPlot(sce.all, reduction = "umap", group.by = "celltype",label = T)
p
ggsave(plot = p, filename='umap_by_celltype.pdf', width = 14, height = 10)
