#################### Step 2.1 Check cell clusters ####################
if(!dir.exists("../2-cell")){
  dir.create("../2-cell")
}
setwd("../2-cell") 

# sce.all <- readRDS("../1-QC/sce.all_best_resolution.Rds")

# For performing differential expression after integration, we switch back to the original data
# DefaultAssay(sce.all) <- "RNA"

# cell types may subject to change according to different projects
# T Cells: CD3D, CD3E
# CD8+ T cells: CD8A, CD8B, CD3D, CD3E
# CD4+ T cells: IL7R, CD3D, CD3E
# B cells: CD19, CD79A, MS4A1 [CD20])
# Plasma cells: IGHG1, MZB1, SDC1, CD79A
# Monocytes and macrophages: CD68, CD163, CD14
# Mast cells: TPSAB1, TPSB2
# NK Cells: NKG7, FGFBP2, FCG3RA, CX3CR1, GZM*
# Photoreceptor cells: RCVRN 
# Fibroblasts: FGF7, MME, CD10
# Endothelial cells: PECAM1, VWF, CD31
# Platelets: PPBP, GP1BB
# Dendritic cells: LILRA4, TPM2
# Monocytes: LYZ, VCAN
# MAIT cells: KLRB1(CD161), CXCR6, CCR6, CCR5, IL7R, IL12R, IL18R, IFNG
# Neutrophils: LTF
# Erythrocytes: HBA1, HBB, HBA2
# Effector and Memory T cells: CCR7-, CD45RA+,CD45RO+
# Epi or tumor: EPCAM, KRT19, PROM1, ALDH1A1, CD24

# Immune (CD45+,PTPRC), Epithelial/cancer (EpCAM+,EPCAM), 
# Stromal (CD10+,MME, Fibo or CD31+,PECAM1,endo)
#

if(T) {
  library(ggplot2) 
  library(tidyverse)
  # Check main cell markers for all major types of cells
  genes_to_check = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A','CD8B','IL7R', # T cells
                     'CD19', 'CD79A', 'MS4A1', # B cells
                     'IGHG1', 'MZB1', 'SDC1', # Plasma cells
                     'CD68', 'CD163', 'CD14', # Mono/Macro cells
                     'TPSAB1' , 'TPSB2',  # Mast cells,
                     'RCVRN','FPR1' , 'ITGAM',
                     'C1QA',  'C1QB',  # Macrophages
                     'S100A9', 'S100A8', 'MMP19','LYZ',# Monocyte
                     'TPM2','LILRA4', # Dendritic cells
                     'LAMP3', 'IDO1','IDO2',## DC3 
                     'CD1E','CD1C', # DC2
                     'KLRB1','NCR1', # NK 
                     'FGF7','MME', 'ACTA2','CD10', # Fibroblasts
                     'PECAM1', 'VWF', # Endothelial cells
                     'HBA1', 'HBB', 'HBA2', # Erythrocytes
                     'EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1') # Epi or tumor cells
  library(stringr)  
  genes_to_check <- convertHumanGeneList(genes_to_check)
  p_all_markers <- DotPlot(sce.all, features = genes_to_check, assay='RNA') + coord_flip()
  
  p_all_markers
  ggsave(plot=p_all_markers, filename="check_all_marker_by_seurat_cluster.pdf",device = "pdf", width = 14, height = 10)
  
  # Check T cell markers
  genes_to_check = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A',
                     'CCR7', 'SELL' , 'TCF7','CXCR6' , 'ITGA1',
                     'CD45RA', 'CD45RO', # Effector and memory T cells (CCR7-)
                     'FOXP3', 'IL2RA',  'CTLA4','PDCD1', # Regulatory T cells
                     'GZMB', 'GZMK','CCL5', # Cytotoxic T cells
                     'IFNG', 'CCL4', 'CCL3' , 
                     'PRF1' , 'NKG7', # NK T cells
                     'KLRB1') # MAIT cells, CXCR6
  genes_to_check <- convertHumanGeneList(genes_to_check)
  
  library(stringr)  
  p <- DotPlot(sce.all, features = genes_to_check, assay='RNA')  + coord_flip()
  p
  ggsave(plot=p, filename="check_Tcells_marker_by_seurat_cluster.pdf", device = "pdf", width = 14, height = 10)
  
  # Check B cell markers
  # mast cells, TPSAB1 and TPSB2 
  # B cell,  CD79A  and MS4A1 (CD20) 
  # naive B cells, such as MS4A1 (CD20), CD19, CD22, TCL1A, and CD83, 
  # plasma B cells, such as CD38, TNFRSF17 (BCMA), and IGHG1/IGHG4
  genes_to_check = c('CD3D','MS4A1','CD79A',
                     'CD19', 'CD22', 'TCL1A',  'CD83', #  naive B cells
                     'CD38','TNFRSF17','IGHG1','IGHG4', # plasma B cells,
                     'TPSAB1' , 'TPSB2',  # Mast cells,
                     'PTPRC' ) 
  genes_to_check <- convertHumanGeneList(genes_to_check)
  p <- DotPlot(sce.all, features = genes_to_check, assay='RNA')  + coord_flip()
  p
  ggsave(plot=p, filename="check_Bcells_marker_by_seurat_cluster.pdf", device = "pdf", width = 14, height = 10)
  
  # Check myeloid cell markers
  genes_to_check = c('CD68', 'CD163', 'CD14',
                     'CD86', 'LAMP3','LILRA4', 'TPM2', #Dendritic cells: LILRA4, TPM2
                     'MRC1','MSR1','ITGAE','ITGAM','ITGAX','SIGLEC7', 
                     'MAF','APOE','FOLR2','RELB','BST2','BATF3')
  genes_to_check <- convertHumanGeneList(genes_to_check)
  
  p <- DotPlot(sce.all, features = unique(genes_to_check), assay='SCT'  )  + coord_flip()
  p
  ggsave(plot=p, filename="check_Myeloid_marker_by_seurat_cluster.pdf", device = "pdf", width = 14, height = 10)
  
  # Check epithelial cell markers
  # epi or tumor (EPCAM, KRT19, PROM1, ALDH1A1, CD24).
  # - alveolar type I cell (AT1; AGER+)
  # - alveolar type II cell (AT2; SFTPA1)
  # - secretory club cell (Club; SCGB1A1+)
  # - basal airway epithelial cells (Basal; KRT17+)
  # - ciliated airway epithelial cells (Ciliated; TPPP3+) 
  
  genes_to_check = c('EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1' ,
                     'AGER','SFTPA1','SCGB1A1','KRT17','TPPP3',
                     'KRT4','KRT14','KRT8','KRT18',
                     'CD3D','PTPRC' ) 
  genes_to_check <- convertHumanGeneList(genes_to_check)
  
  p <- DotPlot(sce.all, features = unique(genes_to_check), assay='RNA')  + coord_flip()
  p
  ggsave(plot=p, filename="check_epi_marker_by_seurat_cluster.pdf", device = "pdf", width = 14, height = 10)
  
  # Check stromal cell markers
  genes_to_check = c('TEK',"PTPRC","EPCAM","PDPN","PECAM1",'PDGFRA','PDGFRB',
                     'CSPG4','GJB2', 'RGS5','ITGA7',
                     'ACTA2','RBP1','CD36', 'ADGRE5','COL11A1','FGF7', 'MME')
  genes_to_check <- convertHumanGeneList(genes_to_check)
  
  p <- DotPlot(sce.all, features = unique(genes_to_check), assay='RNA')  + coord_flip()
  p
  ggsave(plot=p, filename="check_stromal_marker_by_seurat_cluster.pdf", device = "pdf", width = 14, height = 10)
  
  
  p_all_markers
  p_umap=DimPlot(sce.all, reduction = "umap", group.by = "seurat_clusters",label = T, raster=F) 
  
  library(patchwork)
  pau <- p_all_markers+p_umap
  ggsave(plot = pau, filename = 'markers_umap.pdf', device = "pdf", width = 18,height = 12)
  oidP <- DimPlot(sce.all, reduction = "umap",split.by = 'orig.ident', group.by = "seurat_clusters",label = T, raster=F) 
  ggsave(plot=oidP, filename = 'umap_by_orig.ident.pdf', device = "pdf", width = 14, height = 10)
  stdP <- DimPlot(sce.all, reduction = "umap",split.by = 'group', group.by = "seurat_clusters",label = T, raster=F)         
  ggsave(plot = stdP, filename = 'umap_by_group.pdf', device = "pdf", width = 14, height = 10)
}

#################### Step 2.2 Find and visualize top marker genes ####################
# find all markers for each cluster and visualize them.
if(T) {
  # Find marker genes for each group
  sce.markers <- FindAllMarkers(sce.all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  sce.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
  
  # DT::datatable(sce.markers)
  
  pro='cca'
  write.csv(sce.markers,file = paste0(pro,"_sce.all.markers.csv"))
  
  saveRDS(sce.markers, file = paste0(pro,"sce.all.markers.Rds"))
  
  library(dplyr) 
  top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
  DoHeatmap(sce.all,top10$gene,size=3)
  ggsave(filename=paste0(pro,'_sce.markers_heatmap.pdf'),width = 18,height = 24)
  p <- DotPlot(sce.all, features = unique(top10$gene), assay='RNA')  + coord_flip()
  p
  ggsave(plot=p, filename=paste0(pro,'DotPlot_check_top10_markers_by_clusters.pdf'), device="pdf",width = 18,height = 24)
  
  top6 <- sce.markers %>% group_by(cluster) %>% top_n(6, avg_log2FC)
  DoHeatmap(sce.all,top6$gene,size=3)
  ggsave(paste0(pro,'DoHeatmap_check_top6_markers_by_clusters.pdf'),width = 18,height = 12)
  p <- DotPlot(sce.all, features = unique(top6$gene), assay='RNA')  + coord_flip()
  p
  ggsave(plot=p, filename=paste0(pro,'DotPlot_check_top6_markers_by_clusters.pdf'),device="pdf",width = 18,height = 12)
  
  top3 <- sce.markers %>% group_by(cluster) %>% top_n(3, avg_log2FC)
  DoHeatmap(sce.all,top3$gene,size=3)
  ggsave(paste0(pro,'DoHeatmap_check_top3_markers_by_clusters.pdf'),width = 18,height = 12)
  p <- DotPlot(sce.all, features = unique(top3$gene), assay='RNA')  + coord_flip()
  p
  ggsave(plot=p, filename=paste0(pro,'DotPlot_check_top3_markers_by_clusters.pdf'),device="pdf",width = 18,height = 12)
  
  p <- FeaturePlot(sce.all, features = c('MS4A1', 'CD79A','PTPRC','PPBP', 'HBA1', 'HBA2', 'HBB'))
  p
  ggsave(plot=p, filename="check_background_genes.pdf",device="pdf",width = 18,height = 12)
  
}

#################### Step 2.3 Check well-known customized markers ####################
## custom plot for this project
## Feature for major clusters

# 1. Mono Macro: CD14, CD68, LYZ
genes_to_check <- c("CD14", "CD68", "LYZ")
genes_to_check <- convertHumanGeneList(genes_to_check)

p <- FeaturePlot(sce.all, features = genes_to_check, min.cutoff = 1, max.cutoff = 4, raster = F)
p
ggsave(plot=p, filename="MonoMacro_featurePlot.pdf",device="pdf",width = 12,height = 10)


# 2. Granulocyte: FCGR3B, CXCR2, S100A8, S100A9
genes_to_check <- c("FCGR3B", "CXCR2", "S100A8","S100A9")
genes_to_check <- convertHumanGeneList(genes_to_check)

p <- FeaturePlot(sce.all, features = genes_to_check, min.cutoff = 1, max.cutoff = 6, raster = F)
p
ggsave(plot=p, filename="Granulocyte_featurePlot.pdf",device="pdf",width = 12,height = 10)

# 3. T cells: IL7R, IL32, CD3D, CD3G, CD2, CD3E
genes_to_check <- c("IL7R", "IL32", "CD3D","CD3G","CD3E","CD2")
genes_to_check <- convertHumanGeneList(genes_to_check)

p <- FeaturePlot(sce.all, features = genes_to_check, min.cutoff = 1, max.cutoff = 4, raster = F)
p
ggsave(plot=p, filename="Tcells_featurePlot.pdf",device="pdf",width = 12,height = 10)

# 4. B/Plasma: IGHM, MS4A1, CD79A, CD79B, IGHD
genes_to_check <- c("IGHM", "MS4A1", "CD79A","CD79B","IGHD")
genes_to_check <- convertHumanGeneList(genes_to_check)

p <- FeaturePlot(sce.all, features = genes_to_check, min.cutoff = 1, max.cutoff = 5, raster = F)
p
ggsave(plot=p, filename="B_Plasma_cells_featurePlot.pdf",device="pdf",width = 12,height = 10)

# 5. NK cells: GNLY, NKG7, GZMB, KLRD1
genes_to_check <- c("GNLY", "NKG7", "GZMB","KLRD1")
genes_to_check <- convertHumanGeneList(genes_to_check)

p <- FeaturePlot(sce.all, features = genes_to_check, min.cutoff = 1, max.cutoff = 4, raster = F)
p
ggsave(plot=p, filename="NK_cells_featurePlot.pdf",device="pdf",width = 12,height = 10)

# 6. DC cells: FCER1A, CD1C, MRC1
genes_to_check <- c("FCER1A", "CD1C", "MRC1")
genes_to_check <- convertHumanGeneList(genes_to_check)

p <- FeaturePlot(sce.all, features = genes_to_check, min.cutoff = 0.1, max.cutoff = 2, raster = F)
p
ggsave(plot=p, filename="DC_cells_featurePlot.pdf",device="pdf",width = 12,height = 10)

# 7. HSC: GATA1, GATA2
genes_to_check <- c("GATA1", "GATA2")
genes_to_check <- convertHumanGeneList(genes_to_check)

p <- FeaturePlot(sce.all, features = genes_to_check, min.cutoff = 0.1, raster = F)
p
ggsave(plot=p, filename="HSC_cells_featurePlot.pdf",device="pdf",width = 12,height = 10)

##### Check cell markers by cell type
### For bone marrow immune cells
genes_to_check = list(
  T_cells=c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A','CD8B','IL7R'), 
  B_cells=c('CD19', 'CD79A', 'MS4A1'), 
  Plasma_cells=c('IGHG1', 'MZB1', 'SDC1'), 
  Mono_or_Macro=c('CD68', 'CD163', 'CD14'), 
  Mast_cells=c('TPSAB1' , 'TPSB2'),
  Neutrophils=c("MPO","CEACAM8"),
  Eosinophils=c("SIGLEC8"),
  Basophils=c("IL3RA", "FCER1"),
  Macrophages=c('FPR1' , 'ITGAM','C1QA',  'C1QB'),
  Monocyte=c('S100A9', 'S100A8', 'MMP19','LYZ'),
  Dendritic_cells=c('TPM2','LILRA4'), 
  DC3 =c('LAMP3', 'IDO1','IDO2'),
  DC2=c('CD1E','CD1C'),
  Megakaryocytes = c("ITGA2B","GP1BA","GP9","ITGB3","MPL","CXCL4","ITGA2"),
  NK =c('KLRB1','NCR1'),
  Fibroblasts=c('FGF7','MME', 'ACTA2'), 
  Endothelial_cells=c('PECAM1', 'VWF'), 
  HSCs=c("CD34"),
  #Erythrocytes=c('HBA1', 'HBB', 'HBA2'), 
  Epi_or_tumor=c('EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1'))
genes_to_check <- convertHumanGeneList(genes_to_check)

p_all_markers <- DotPlot(sce.all, features = genes_to_check, assay='RNA') +
  theme(axis.text.x=element_text(angle=45,hjust=1.2,vjust = 1.1,size = 8),strip.text = element_text(size = 8) )

p_all_markers

ggsave("check_all_cell_markers_from_BM.pdf", plot = p_all_markers, width = 14, height = 10)

#### For T cell subcluster
genes_to_check = list(
  Mono_Macro=c('CD68', 'CD163', 'CD14'), 
  T_cells=c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A','CD8B'), 
  gdT_cells=c('TRDC','TRGC2','CMC1', 'ITGAD'), 
  regT = c("FOXP3","CTLA4","IL2RA"),
  Tex = c('PDCD1','LAG3', 'HAVCR2','TIGIT','TOX'),
  naive_T=c("SELL","CCR7","IL7R","LTB"), 
  Tc_cells = c('CX3CR1', 'LGALS1', 'FGFBP2', 'GNLY', 'GZMB', 'NKG7'),
  Tem_cells = c("CD52","TNFRSF4","AQP3","TIMP1",'KLRG1','ITGA1','CXCR3'),
  NKT =c("TRAV1-1","TRBV25-1"),
  NK_cells =c('KLRB1','NCR1'),
  Macrophages=c('FPR1' , 'ITGAM','C1QA',  'C1QB'),
  Monocyte=c('S100A9', 'S100A8', 'MMP19','LYZ'),
  Fibroblasts=c('FGF7','MME', 'ACTA2'), 
  Endothelial_cells=c('PECAM1', 'VWF'), 
  Erythrocytes=c('HBA1', 'HBB', 'HBA2'), 
  Epi_or_tumor=c('EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1'))

p_all_markers <- DotPlot(sce, features = genes_to_check, assay='RNA') +
  theme(axis.text.x=element_text(angle=45,hjust=1.2,vjust = 1.1,size = 8),strip.text = element_text(size = 8) )

p_all_markers

ggsave("check_all_cell_markers_for_T.pdf", plot = p_all_markers, width = 14, height = 10)


#### For mouse T cell subcluster
genes_to_check = list(
  Mono_Macro= c("Cd68","Cd163", "Cd14"), # c('CD68', 'CD163', 'CD14'), 
  T_cells= c("Il7r","Cd8a","Cd8b1", "Cd4", "Cd3d","Cd3e"), # c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A','CD8B')
  gdT_cells=c("Trdc", "Cmc1", "Itgad"), #c('TRDC','TRGC2','CMC1', 'ITGAD'), 
  regT = c("Foxp3", "Ctla4", "Il2ra"), #c("FOXP3","CTLA4","IL2RA"),
  Tex = c("Tigit", "Lag3", "Tox", "Havcr2","Pdcd1"),  # c('PDCD1','LAG3', 'HAVCR2','TIGIT','TOX'),
  naive_T= c("Ccr7", "Ltb", "Sell"), #c("SELL","CCR7","IL7R","LTB"), 
  Tc_cells = c("Nkg7","Gzme","Gzmd","Gzmg","Gzmn","Gzmf","Gzmc","Gzmb","Lgals1","Cx3cr1"), # c('CX3CR1', 'LGALS1', 'FGFBP2', 'GNLY', 'GZMB', 'NKG7'),
  Tem_cells = c("Klrg1","Timp1","Cd52","Cxcr3","Tnfrsf4","Itga1","Aqp3" ), #c("CD52","TNFRSF4","AQP3","TIMP1",'KLRG1','ITGA1','CXCR3'),
  NKT = c("Trav1"), # c("TRAV1-1","TRBV25-1"),
  NK_cells = c("Ncr1","Klrb1a","Klrb1","Klrb1c","Gm44511","Klrb1b", "Klrb1f"), # c('KLRB1','NCR1'),
  Macrophages= c("C1qb","Fpr1","C1qa","Itgam","Gm49368"), # c('FPR1' , 'ITGAM','C1QA',  'C1QB'),
  Monocyte= c("9530003J23Rik","Gm5849","Mmp19","S100a8"), # c('S100A9', 'S100A8', 'MMP19','LYZ'),
  Fibroblasts= c("Acta2", "Mme","Fgf7" ), # c('FGF7','MME', 'ACTA2'), 
  Endothelial_cells= c("Vwf","Pecam1"), # c('PECAM1', 'VWF'), 
  Erythrocytes= c("Hba-a1", "Hba-a2"), # c('HBA1', 'HBB', 'HBA2'), 
  Epi_or_tumor= c("Epcam","Aldh1a1", "Aldh1a7", "Prom1", "Krt19" )) # c('EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1'))

p_all_markers <- DotPlot(sce.all, features = genes_to_check, assay='RNA') +
  theme(axis.text.x=element_text(angle=45,hjust=1.2,vjust = 1.1,size = 8),strip.text = element_text(size = 8) )

p_all_markers

ggsave("check_all_cell_markers_for_T.pdf", plot = p_all_markers, width = 14, height = 10)


DotPlot(sce.all, features = c("Il7r", "Lef1", "Tcf7", "Cxcr6", "Satb1", "Tet2", "Sell", "Ccr7","Tox","Il12b"), group.by = "seurat_clusters") + coord_flip()

DotPlot(sce.all, features = c("Il12b", "Cxcr4", "Il6", "Foxp3"), group.by = "seurat_clusters") + coord_flip()




### For pancreatic cells
genes_to_check = list(
  Beta=c("Ins1","Ins2"),
  Alpha=c("Gcg"),
  Delta=c("Sst"),
  Acinar=c("Cpa1"),
  Ductal=c("Krt17"),
  PP=c("Ppy"),
  ProlifEndo=c("Gmnn","Sst","Gcg","Ppy"),
  Mesenchymal=c("Col1a"),
  Macrophages=c("Ccr5"),
  ECs=c("Pecam1","Cd32"),
  Others=c("Pparg","Smad7","Wnt5a"))

p_all_markers <- DotPlot(sce.all, features = genes_to_check, assay='RNA') +
  theme(axis.text.x=element_text(angle=45,hjust=1.2,vjust = 1.1,size = 8),strip.text = element_text(size = 8) )

p_all_markers

ggsave("check_all_cell_markers_from_pancrease.pdf", plot = p_all_markers, width = 14, height = 10)

#### For Moncyte and Macrophage subcluster
genes_to_check = list(
  Mono_Macro = c('CD68', 'CD163', 'CD14', "FCN1"), 
  CD14_Mono = c("S100A8","S100A9","VCAN"),
  CD16_Mono = c("FCGR3A","LILRB2","LST1"),
  LYVE1_Macro = c("LYVE1","C1QC","PLTP","SEPP1"),
  NLRP3_Macro = c('NLRP3','IL1B','CXCL2','EREG'),
  Alveolar_Marco = c('PPARG','MARCO','MRC1','MSR'),
  ISG15_Macro = c('ISG15','CXCL10','IFITM3','GBP1'),
  C1QC_Macro = c('C1QA','C1QB','APOE'),
  SPP1_Macro = c('SPP1','VEGFA','GPNMB','FN1'),
  T_cells = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A','CD8B'), 
  gdT_cells = c('TRDC','TRGC2','CMC1', 'ITGAD'), 
  regT = c("FOXP3","CTLA4","IL2RA"),
  Tex = c('PDCD1','LAG3', 'HAVCR2','TIGIT','TOX'),
  naive_T = c("SELL","CCR7","IL7R","LTB"), 
  Tc_cells = c('CX3CR1', 'LGALS1', 'FGFBP2', 'GNLY', 'GZMB', 'NKG7'),
  Tem_cells = c("CD52","TNFRSF4","AQP3","TIMP1",'KLRG1','ITGA1','CXCR3'),
  NKT = c("TRAV1-1","TRBV25-1"),
  NK_cells = c('KLRB1','NCR1'),
  Macrophages = c('FPR1' , 'ITGAM'),
  Monocyte = c('MMP19','LYZ'),
  Fibroblasts = c('FGF7','MME', 'ACTA2'), 
  Endothelial_cells = c('PECAM1', 'VWF'), 
  Erythrocytes = c('HBA1', 'HBB', 'HBA2'), 
  Epi_or_tumor = c('EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1'))

p_all_markers <- DotPlot(sce, features = genes_to_check, assay='RNA') +
  theme(axis.text.x=element_text(angle=45,hjust=1.2,vjust = 1.1,size = 8),strip.text = element_text(size = 8) )

p_all_markers

ggsave("check_all_cell_markers_for_MonoMacro.pdf", plot = p_all_markers, width = 14, height = 10)