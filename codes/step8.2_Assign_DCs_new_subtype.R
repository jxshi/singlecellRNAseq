############### Step 8.2 Assign DC cell newsubtype.

setwd("/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype/DCs/")

sce.sub <- readRDS("DCs.sce.sub.harmony.Rds")
sce.sub$group <- factor(sce.sub$group, levels = c("PBS", "DD_mGE", "DD_mIL12", "DD_mGE12", "TD_mGE12"))
sce.all$group <- factor(sce.all$group, levels = c("PBS", "DD_mGE", "DD_mIL12", "DD_mGE12", "TD_mGE12"))
## Dendritic cell subtype marker panel
dc_markers <- list(
  DC_general        = c("Itgax", "H2-Ab1"),
  moDC_cDC2         = c("Ly6c2", "Fcgr1", "Mki67"),
  Th2_polarizing_cDC2 = c("Ccl17", "Mgl2", "Mrc1"),
  cDC1              = c("Xcr1", "Clec9a", "Batf3", "Irf8", "Itgae", "Cd8a"),
  cDC2              = c("Itgam", "Sirpa", "Irf4", "Fcgr1"),
  migratory_cDC2    = c("Ccr7"),
  pDC               = c("Bst2", "Siglech", "Il3ra", "Tlr7", "Tlr9"), # Itgax low
  T                 = c("Cd3d","Cd3e", "Cd4", "Cd8a","Cd8b1") 
)
DotPlot(sce.sub, features = unique(unlist(dc_markers)), group.by = "newsubcelltype") + RotatedAxis()
DotPlot(sce.sub, features = unique(unlist(dc_markers)), group.by = "seurat_clusters") + RotatedAxis()
Idents(sce.sub) <- sce.sub$seurat_clusters
if(T) {
  table(sce.sub@active.ident)
  
  newnewsubcelltype=data.frame(ClusterID=0:9, newsubcelltype='na')
  
  newsubcelltype[newsubcelltype$ClusterID %in% c( 0),2]='moDC'
  newsubcelltype[newsubcelltype$ClusterID %in% c( 1),2]='cDC2'          # chatGPT: cDC2; Grok: cDC2; Gemini: moDC or TAM; 
  newsubcelltype[newsubcelltype$ClusterID %in% c( 2),2]='cDC1'          # chatGPT: cDC1; Grok: cDC1; Gemini: cDC1;
  newsubcelltype[newsubcelltype$ClusterID %in% c( 3),2]='pDC'           # chatGPT: pDC; Grok: pDC; Gemini: pDC;
  newsubcelltype[newsubcelltype$ClusterID %in% c( 4),2]='migratoryDC'   # chatGPT: Fscn1⁺ Ccr7⁺ Ccl17⁺/Ccl22⁺ migratory DC; Grok: migratoryDC; Gemini: migratoryDC;
  newsubcelltype[newsubcelltype$ClusterID %in% c( 5),2]='moDC'           # chatGPT: Ly6C⁺ CCR2⁺ pre-migratory mo-DC; Grok: pre-DC  ;Gemini: cDC2
  newsubcelltype[newsubcelltype$ClusterID %in% c( 6),2]='CD8 T cells'   # Proliferating T, Cd3e. # 同时表达DC和T的marker基因，是个混群，去除。
  newsubcelltype[newsubcelltype$ClusterID %in% c( 7),2]='migratoryDC'          # chatGPT: CCR7⁺ Fscn1⁺ migratory cDC
  newsubcelltype[newsubcelltype$ClusterID %in% c( 8),2]='Th1'           # chatGPT: Activated Th1-like T
  newsubcelltype[newsubcelltype$ClusterID %in% c( 9),2]='pDC'           # 

  
  
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
  
  maincelltype=data.frame(ClusterID=0:9, maincelltype='na')
  
  maincelltype[maincelltype$ClusterID %in% c( 0),2]='DCs'
  maincelltype[maincelltype$ClusterID %in% c( 1),2]='DCs'
  maincelltype[maincelltype$ClusterID %in% c( 2),2]='DCs'
  maincelltype[maincelltype$ClusterID %in% c( 3),2]='DCs'
  maincelltype[maincelltype$ClusterID %in% c( 4),2]='DCs'
  maincelltype[maincelltype$ClusterID %in% c( 5),2]='DCs'
  maincelltype[maincelltype$ClusterID %in% c( 6),2]='Unknown'
  maincelltype[maincelltype$ClusterID %in% c( 7),2]='DCs'
  maincelltype[maincelltype$ClusterID %in% c( 8),2]='Unknown'
  maincelltype[maincelltype$ClusterID %in% c( 9),2]='DCs'

  
  
  
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
sce.dc <- sce.sub

table(sce.dc$group,sce.dc$maincelltype)

saveRDS(sce.dc, "sce.DCs_w_correction.Rds")

#################################### Plot Figures for DCs ###############################################
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



setwd("/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype/graphs/")
dir.create("DCsFigures")

setwd("DCsFigures/")
# sce.dc <- readRDS("~/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype/DCs/sce.DCs_w_correction.Rds")
# sce.new <- sce.dc
Idents(sce.dc) <- sce.dc$newsubcelltype
sce.new <- subset(sce.dc, idents = setdiff(levels(Idents(sce.dc)), c("CD8 T cells", "Th1")))

Idents(sce.new) <- sce.new$newsubcelltype
sce.new$subcelltype <- sce.new$newsubcelltype


p <- DimPlot(sce.new, group.by = "maincelltype", raster = F, label = T, cols = combined_colors)
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



