library(ggsci)
source("~/software/functions/groupbarcharts.R")
source("~/software/functions/groupbarchartssubcluster.R")
source("~/software/functions/PropPlot.R")
source("~/software/functions/plot_integrated_cluster.R")
default_colors <- pal_nejm(palette = c("default"), alpha=1)(8)
extra_colors <- pal_aaas()(8)
combined_colors <- c(default_colors, extra_colors)
combined_colors

combined_colors <- c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#008B45FF", "#5F559BFF", "#6F99ADFF","#A20056FF",  "#FF5733",
                     "#EE4C97FF", "#3B4992FF", "#008280FF","#631879FF", "#C70039","#FFDC91FF")

combined_colors <- c("#BC3C29", "#0072B5", "#E18727", "#008B45", "#5F559B", "#6F99AD", "#A20056", "#EE4C97", "#3B4992", "#008280", "#631879", "#C70039", "#FFDC91")

combined_colors <- c(
  "#1A1B4B",  # deep navy
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
  "#A62071",  # magenta
  "#C0261D",  # red
  "#D2A21B",  # mustard yellow
  "#F5C34D",  # gold
  "#F1A800",  # bright yellow
  "#DDA873",  # peach-orange
  "#3C8A32"   # medium green
)

combined_colors <- c("#A9B4D7","#9BC9D9", "#F1A800", "#5EACA7", "#BE84B8", "#3C8A32","#D8A5CE")

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

cnames <- readRDS("celltype_name.Rds")
fnames <- readRDS("celltype_file_folder.Rds")
fnames
cnames

setwd("/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype/graphs/")
dir.create("mainfigures")

setwd("mainfigures/")
sce.all <- readRDS("~/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype/sce.all_w_corrected_new_subcelltype_main_celltype_wo_unknown.Rds")
table(sce.all$maincelltype)

sce.new <- subset(sce.all, idents = setdiff(levels(Idents(sce.all)), c("Epithelial cells","Fibroblasts")))

Idents(sce.new) <- sce.new$maincelltype
sce.new$subcelltype <- sce.new$newsubcelltype

sce.new$group <- factor(sce.new$group, levels = c("PBS", "DD_mGE", "DD_mIL12", "DD_mGE12", "TD_mGE12"))

p <- DimPlot(sce.new, group.by = "maincelltype", raster = F, label = T, cols = maincelltype_colors, pt.size =0.0002, shuffle = T, repel = T) 
p  
ggsave("main_celltype_UMAP.pdf", plot = p, device = cairo_pdf, width = 12, height = 12)

p <- DimPlot(sce.new, group.by = "group", raster = F, label = F, cols = group_colors, split.by = "group")
# p <- LabelClusters(plot = p, id = "group", repel = TRUE)
p
ggsave("UMAP_by_group.pdf", plot = p, device = cairo_pdf, width = 12, height = 12)


p <- groupbarcharts(sce.new, "group") 
p
ggsave("compare_maincelltype_percentage_of_each_group.pdf", plot = p, device = cairo_pdf, width = 8, height = 8)
tmp <- table(sce.new$group,sce.new$maincelltype)
write.table(tmp, file = "maincelltype_per_group.tsv", quote = F, col.names = NA, sep = "\t")

p <- groupbarchartssubcluster(sce.new, "group")
p
ggsave("compare_percentage_of_celsl_in_each_group.pdf", plot = p, device = cairo_pdf, width = 24, height = 8)

p <-  plot_integrated_clusters(sce.new, cluster_col = "maincelltype", group_col = "group")
p
ggsave("compare_percentage_of_maincelltype.pdf", plot = p, device = cairo_pdf, width = 8, height = 8)

Idents(sce.new) <- sce.new$maincelltype
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



genes_to_check <- c("Ptprc",'Cd3d','Cd3e','Cd3g','Flt3','H2-Eb1','H2-Aa','A2-Ab1','Cd74')


genes_to_check = list(
  T_cells = c("Ptprc",'Cd3d','Cd3e','Cd3g'), 
  DCs = c('H2-Eb1','H2-Aa','H2-Ab1','Cd74','Flt3',"Clec12a","Clec9a"),
  MonoMacro = c('Cd68','Mrc1','Csf1r','Mafb','C1qa','Itgam'),
  NK = c('Klrb1','Ncr1','Klrc1','Klra8','Klra3','Prf1',"Xcr1"),
  Granulocytes = c('S100a8','S100a9','Ly6g','Ly6e'),
  B_cells = c('Cd79a','Cd79b','Igkc','Ighm','Ms4a1','Fcmr'),
  Mast_cells =c('Tpsab1','Mcpt1','Mcpt2','Tpsb2','Cma1','Cpa3','Fcer1a'))

genes_to_check = list(
  T_cells = c("Ptprc",'Cd3d','Cd3e','Cd3g'), 
  DCs = c('H2-Eb1','H2-Aa','H2-Ab1','Cd74','Flt3'),
  MonoMacro = c('Cd68','Csf1r','Mafb','Itgam','C1qa','Mrc1'),
  NK = c('Ncr1','Klrc1','Prf1','Klra3'),
  Granulocytes = c('S100a8','S100a9','Ly6e'),
  B_cells = c('Cd79a','Igkc','Ighm','Ms4a1','Fcmr'),
  Mast_cells =c('Tpsab1','Cma1','Tpsb2','Cpa3','Mcpt1','Fcer1a'))

p_all_markers <- DotPlot(sce.new, features = genes_to_check, assay='RNA') +
  theme(axis.text.x=element_text(angle=45,hjust=1.2,vjust = 1.1,size = 8),strip.text = element_text(size = 8) )

p_all_markers

source("~/software/functions/PropPlot.R")

p <- PropPlot(sce.new, "group", maincelltype_colors)
p
ggsave("compare_percentage_of_maincelltype_for_each_group.pdf", plot = p, device = cairo_pdf, width = 8, height = 8)


p <- DimPlot(sce.new, group.by = "celltype", raster = F, label = T, cols = maincelltype_colors, pt.size =0.002, shuffle = T, repel = T) + 
  theme(legend.position = "none")
coords <- Embeddings(sce.new, "umap")[, 1:2] %>% as.data.frame()
coords$ident <- Idents(sce.new)

# Add one contour per identity
for (celltype in levels(coords$ident)) {
  p <- p + geom_density_2d(data = subset(coords, ident == celltype),
                           aes(x = UMAP_1, y = UMAP_2),
                           color = scales::hue_pal()(length(levels(coords$ident)))[which(levels(coords$ident) == celltype)],
                           bins = 8, size = 1)
}

p


p <- DimPlot(sce.new, group.by = "celltype", raster = F, label = T, cols = maincelltype_colors, pt.size =0.002, shuffle = T, repel = T)

df <- data.frame(Embeddings(sce.new, "umap"), celltype = sce.new$celltype)

# Add contours per cell type
for (ct in unique(df$celltype)) {
  p <- p + geom_density_2d(
    data = df[df$celltype == ct, ],
    aes(
      x = UMAP_1,
      y = UMAP_2,
      group = celltype,
      color = celltype        # <-- must map to celltype, not colors vector
    ),
    size = 0.4,
    bins = 6
  )
}

p
ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = celltype), size = 0.2, alpha = 0.6) +
  geom_density_2d(aes(color = celltype, group = celltype),
                  size = 0.4, bins = 6) +
  scale_color_manual(values = maincelltype_colors) +
  theme_classic()
