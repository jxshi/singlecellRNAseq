############### Step 8.3 Assign MonoMacro cell newsubtype.
setwd("/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype/MonoMacro/")

sce.sub <- readRDS("MonoMacro.sce.sub.harmony.Rds")

genes_to_check <- c("Cd14","Cd16","Ly6c2","Ly6c1",'Cd68',"Adgre1","Csf1r",'Mertk',"Axl","CD163","Clec7a","Clec10a",'Cd209')
DotPlot(sce.sub, features = genes_to_check)

sce.sub <- sce.monomacro
Idents(sce.sub) <- sce.sub$seurat_clusters

m0_genes <- c("Csf1r","Cd68","Aif1","Lyz2","Adgre1","Mertk","Tyrobp",
              "C1qa","C1qb","C1qc","Lgals3","Cd36","Ccl5","H2-Aa","H2-Ab1")

m1_genes <- c("Nos2","Il1b","Tnf","Il6","Cxcl9","Cxcl10","Stat1","Irf5",
              "Cd86","Cd80","Gbp2","Gbp5", # "H2-Aa","H2-Ab1",
              "Socs3","Ifit1","Ifit3","Isg15","Rsad2","Irf7")

m2_genes <- c("Mrc1","Cd163","Arg1","Chil3","Retnla","Clec10a","Il10",
              "Tgfb1","Mgl2","Ccl17","Ccl22","Vegfa","Mmp12","Spp1",
              "Apoe","Pdgfb","Timp2","Fn1","Anxa1","Stab1","Vsig4",
              "Klf4", "Stat6","Irf4","Cd63","Il10ra","Il4ra","Ly6c1","Ly6c2","Slfn4",'Cd274','Acp5','Ctsk')

DotPlot(sce.sub, features = list(M0 = m0_genes, 
                             M1 = m1_genes, 
                             M2 = m2_genes)) + RotatedAxis()

sce.sub <- AddModuleScore(sce.sub, features = list(M0 = m0_genes), name = "M0_Score")
sce.sub <- AddModuleScore(sce.sub, features = list(M1 = m1_genes), name = "M1_Score")
sce.sub <- AddModuleScore(sce.sub, features = list(M2 = m2_genes), name = "M2_Score")

FeaturePlot(sce.sub, features = c("M0_Score1","M1_Score1","M2_Score1"))
VlnPlot(sce.sub, features = c("M0_Score1","M1_Score1","M2_Score1"), group.by = "seurat_clusters")





if(T) {
  table(sce.sub@active.ident)
  
  newsubcelltype=data.frame(ClusterID=0:9, newsubcelltype='na')
  
  newsubcelltype[newsubcelltype$ClusterID %in% c( 0),2]='Macrophages0'  # highly activated, inflammatory macrophages (M1-like TAMs)
  newsubcelltype[newsubcelltype$ClusterID %in% c( 1),2]='Macrophages1'  # Ly6C⁺ monocyte-derived TAM
  newsubcelltype[newsubcelltype$ClusterID %in% c( 2),2]='Macrophages2'  # IFN/ISG-high inflammatory TAM; Gemini: M1-like TAMs
  newsubcelltype[newsubcelltype$ClusterID %in% c( 3),2]='Macrophages3'  # MHC-II+ Cxcl9⁺ TAM
  newsubcelltype[newsubcelltype$ClusterID %in% c( 4),2]='Macrophages4'  # OXPHOS-high Macrophage 
  newsubcelltype[newsubcelltype$ClusterID %in% c( 5),2]='Macrophages5'  # C1q⁺ M2
  newsubcelltype[newsubcelltype$ClusterID %in% c( 6),2]='Macrophages6'  # Slfn4⁺ Pd-L1⁺ M2
  newsubcelltype[newsubcelltype$ClusterID %in% c( 7),2]='Macrophages7'  # Mrc1+ Arg1+ M2. tissue-remodeling / pro-angiogenic TAMs
  newsubcelltype[newsubcelltype$ClusterID %in% c( 8),2]='Macrophages8'  # Arg1⁺ Spp1⁺ hypoxia/glycolysis-high M2; tissue-remodeling / pro-angiogenic TAMs
  newsubcelltype[newsubcelltype$ClusterID %in% c( 9),2]='Macrophages9'  # osteoclast-like macrophage
  
  
  
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
  Idents(sce.sub) <- sce.sub$seurat_clusters
  
  maincelltype=data.frame(ClusterID=0:9, maincelltype='na')
  
  maincelltype[maincelltype$ClusterID %in% c( 0),2]='Macrophages'
  maincelltype[maincelltype$ClusterID %in% c( 1),2]='Macrophages'
  maincelltype[maincelltype$ClusterID %in% c( 2),2]='Macrophages'
  maincelltype[maincelltype$ClusterID %in% c( 3),2]='Macrophages'
  maincelltype[maincelltype$ClusterID %in% c( 4),2]='Macrophages'
  maincelltype[maincelltype$ClusterID %in% c( 5),2]='Macrophages'
  maincelltype[maincelltype$ClusterID %in% c( 6),2]='Macrophages'
  maincelltype[maincelltype$ClusterID %in% c( 7),2]='Macrophages'
  maincelltype[maincelltype$ClusterID %in% c( 8),2]='Macrophages'
  maincelltype[maincelltype$ClusterID %in% c( 9),2]='Macrophages'
  
  
  
  
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
dir.create("MonoMacroFigures")

setwd("MonoMacroFigures/")
# sce.dc <- readRDS("~/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype/DCs/sce.DCs_w_correction.Rds")
sce.new <- sce.monomacro

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
