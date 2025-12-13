setwd("/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype")

#################### Step 3.4 Split data according to modified cell types ####################
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)

sce.all <- readRDS("sce.new_w_corrected_main_celltype.Rds")

#拆分为 多个 seurat子对象
sce.all.list <- SplitObject(sce.all, split.by = "maincelltype")
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


genes_to_check <- c("Cd14", "Cd16","Marco")
FeaturePlot(sce.all, genes_to_check)
