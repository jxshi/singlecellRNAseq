setwd("/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/3-celltype/Bcells")
sce.sub <- readRDS("Bcells.sce.sub.harmony.Rds")
if(T) {
  table(sce.sub@active.ident)
  
  subcelltype=data.frame(ClusterID=0:15, subcelltype='na')
  
  subcelltype[subcelltype$ClusterID %in% c( 0),2]='germinal center like B'
  subcelltype[subcelltype$ClusterID %in% c( 1),2]='memory like B'
  subcelltype[subcelltype$ClusterID %in% c( 2),2]='naive B'
  subcelltype[subcelltype$ClusterID %in% c( 3),2]='Unknown'
  subcelltype[subcelltype$ClusterID %in% c( 4),2]='APC-B'
  subcelltype[subcelltype$ClusterID %in% c( 5),2]='pre plasmablast B'
  subcelltype[subcelltype$ClusterID %in% c( 6),2]='proliferating B'
  subcelltype[subcelltype$ClusterID %in% c( 7),2]='Unknown'
  subcelltype[subcelltype$ClusterID %in% c( 8),2]='germinal center like B'
  subcelltype[subcelltype$ClusterID %in% c( 9),2]='Unknown'
  subcelltype[subcelltype$ClusterID %in% c(10),2]='Unknown'
  subcelltype[subcelltype$ClusterID %in% c(11),2]='transitional B'
  subcelltype[subcelltype$ClusterID %in% c(12),2]='Unknown'
  subcelltype[subcelltype$ClusterID %in% c(13),2]='Unknown'
  subcelltype[subcelltype$ClusterID %in% c(14),2]='APC-B'
  subcelltype[subcelltype$ClusterID %in% c(15),2]='germinal center like B'
  
  
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
  
  maincelltype[maincelltype$ClusterID %in% c( 0),2]='B cells'
  maincelltype[maincelltype$ClusterID %in% c( 1),2]='B cells'
  maincelltype[maincelltype$ClusterID %in% c( 2),2]='B cells'
  maincelltype[maincelltype$ClusterID %in% c( 3),2]='Unknown'
  maincelltype[maincelltype$ClusterID %in% c( 4),2]='B cells'
  maincelltype[maincelltype$ClusterID %in% c( 5),2]='B cells'
  maincelltype[maincelltype$ClusterID %in% c( 6),2]='B cells'
  maincelltype[maincelltype$ClusterID %in% c( 7),2]='Unknown'
  maincelltype[maincelltype$ClusterID %in% c( 8),2]='B cells'
  maincelltype[maincelltype$ClusterID %in% c( 9),2]='Unknown'
  maincelltype[maincelltype$ClusterID %in% c(10),2]='Unknown'
  maincelltype[maincelltype$ClusterID %in% c(11),2]='B cells'
  maincelltype[maincelltype$ClusterID %in% c(12),2]='Unknown'
  maincelltype[maincelltype$ClusterID %in% c(13),2]='Unknown'
  maincelltype[maincelltype$ClusterID %in% c(14),2]='B cells'
  maincelltype[maincelltype$ClusterID %in% c(15),2]='B cells'
  
  
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
sce.b <- sce.sub

table(sce.b$group,sce.b$maincelltype)

saveRDS(sce.b, "sce.Bcells_w_correction.Rds")
