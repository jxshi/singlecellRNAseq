### ---------------
###
### Author: Jianxiang Shi
### Date: 2022-11-05 22:28:00
### BGI College & Henan Institute of Medical and Pharmaceutical Sciences, Zhengzhou University
### Update Log: 2022-02-08  First version
### Update Log: 2022-11-05  Second version
###
### ---------------

################################# Put your project information here #################################
## Project Name: Put a descriptive project NAME here.
## Collabortor: 
## Original Data Location:
## ~/data/singlecell/bgi/wangpengju/xuanyujing/results: lw5, lwcd45. MM CD45+  cells, relapse after CAR-T treatment.

#####################################################################################################

rm(list = ls())
options(stringsAsFactors = F)

# Set up the working directory here. ****** IMPORTANT ******
setwd("~/data/singlecell/bgi/wangpengju/xuanyujing/results")

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

## Parallellization in Seurat with future
## https://satijalab.org/seurat/articles/future_vignette.html
# plan("multicore", workers = 10)

# Ref1: https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1
## Ref2: 使用参考样本进行的多中心基准测试scRNA-seq技术的研究
## Ref2 : https://mp.weixin.qq.com/s/TnaXVYlUMVXLC4qIZyZ0OA
## https://swaruplab.bio.uci.edu/tutorial/integration/integration_tutorial.html

#################### Step 1. Load data ####################
# If data has been loaded previously, load previously saved Rds ojbect directly.
if(file.exists("sceList.raw.Rds")) {
  sceList <- readRDS("sceList.raw.Rds")
  
  samples = list.files('../datamatrix/')
  samples
  length(samples)
}

# Otherwise, read data from data matrix for 10X genomics data
if(!file.exists("sceList.raw.Rds")) {
  
  samples = list.files('../datamatrix/')
  samples
  length(samples)
  
  ## set default filtering parameters to filter empty cells
  sceList = lapply(samples,function(pro){ 
    sce <- CreateSeuratObject(counts = Read10X(paste0("../datamatrix/", pro,"/")), 
                              min.cells = 3, 
                              min.features = 200,
                              project = pro)
    return(sce)
  })
  
  names(sceList)
  
  group <- sub("^([^_]+_[^_]+)_.*", "\\1", samples)
  group
  
  group[10:12] <- "PBS"
  
  table(group)
  
  for(i in seq(length(sceList))) {
    sceList[[i]]@project.name <- samples[i]
    sceList[[i]][['batch']] <- samples[i]
    sceList[[i]][['sample']] <- samples[i]
    sceList[[i]][['group']] <- group[i]
    levels(sceList[[i]]@active.ident) <- samples[i]
  }
  
  saveRDS(sceList,file="sceList.raw.Rds")
}

#################### Step 2. Quality Control ####################
#dir.create("./1-QC")
setwd("./1-QC")

# Visualization of mitochondrial, ribosomal, homoglobin percent before and after filtering
## Ref: 参考耗大爷的单细胞RNA-Seq数据分析Manual
## SCTransform in included in this step
## Ref: https://satijalab.org/seurat/articles/sctransform_vignette.html
if(file.exists("sceList.Rds")) {
  sceList <- readRDS("sceList.Rds")
}
# Calculate mitochondrial, ribosomal, homoglobin percent
if(!file.exists("sceList.Rds")) {
  
  set.seed(010101)
  
  # Cell cycle genes for mouse research
  # s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
  # g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)
  
  # Read from local directory.
  s.genes <- readRDS("~/software/functions/m.s.genes.Rds")
  g2m.genes <- readRDS("~/software/functions/m.g2m.genes.Rds")
  
  sceList <- lapply(X = sceList, FUN = function(x) {
    x <- PercentageFeatureSet(x, pattern = "^[Mm]t-", col.name = "percent.mito") %>%
      PercentageFeatureSet(pattern = "^Rp[sl]", col.name = "percent.ribo") %>%
      PercentageFeatureSet(pattern = "^Hb[^(p)]", col.name = "percent.hb") %>%
      NormalizeData() %>%
      CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) %>%
      SCTransform(method = "glmGamPoi", do.correct.umi = TRUE, vars.to.regress = c("nCount_RNA","S.Score", "G2M.Score"),
                  do.scale = TRUE, do.center = TRUE, verbose = FALSE)
  })
  
  for (i in seq(length(sceList))){
    Idents(sceList[[i]]) <- sceList[[i]]$sample
  }
  
  sceList.old <- sceList
  
  for (i in seq(length(sceList))) {
    sceList.old[[i]][['filtered']] <- 'non-filtered'
    
    sceList[[i]] <- subset(sceList[[i]], subset=nFeature_RNA > 200 & percent.mito < 20 & percent.hb < 10)
    sceList[[i]][['filtered']] <- 'filtered'
  }
  
  cl <- makeCluster(getOption('cl.cores', length(samples)))
  sceList.merged <- list(merge(sceList.old[[1]], sceList[[1]]),
                         merge(sceList.old[[2]], sceList[[2]]),
                         merge(sceList.old[[3]], sceList[[3]]),
                         merge(sceList.old[[4]], sceList[[4]]),
                         merge(sceList.old[[5]], sceList[[5]]),
                         merge(sceList.old[[6]], sceList[[6]]),
                         merge(sceList.old[[7]], sceList[[7]]),
                         merge(sceList.old[[8]], sceList[[8]]),
                         merge(sceList.old[[9]], sceList[[9]]),
                         merge(sceList.old[[10]], sceList[[10]]),
                         merge(sceList.old[[11]], sceList[[11]]),
                         merge(sceList.old[[12]], sceList[[12]])
  )
  
  pList1 <- parLapply(cl, sceList.merged, VlnPlot, features=c('nFeature_RNA', 'nCount_RNA', 'percent.mito', 'percent.hb'), pt.size=0.1, split.by='filtered', split.plot=TRUE, combine=FALSE)
  pList2 <- parLapply(cl, sceList.merged, FeatureScatter, feature1='nCount_RNA', feature2='percent.mito', pt.size=0.1, group.by='filtered')
  pList3 <- parLapply(cl, sceList.merged, FeatureScatter, feature1='nCount_RNA', feature2='nFeature_RNA', pt.size=0.1, group.by='filtered')
  pList4 <- parLapply(cl, sceList.merged, FeatureScatter, feature1='nCount_RNA', feature2='percent.ribo', pt.size=0.1, group.by='filtered')
  pList5 <- parLapply(cl, sceList.merged, FeatureScatter, feature1='nCount_RNA', feature2='percent.hb', pt.size=0.1, group.by='filtered')
  
  subP1 <- (pList1[[1]][[1]] + FontSize(x.title=0)) + (pList1[[2]][[1]] + FontSize(x.title=0))
  subP2 <- (pList1[[1]][[2]] + FontSize(x.title=0)) + (pList1[[2]][[2]] + FontSize(x.title=0))
  subP3 <- (pList1[[1]][[3]] + FontSize(x.title=0)) + (pList1[[2]][[3]] + FontSize(x.title=0))
  
  subP4 <- (pList2[[1]] + RotatedAxis()) + (pList2[[2]] + RotatedAxis())
  subP5 <- (pList3[[1]] + RotatedAxis()) + (pList3[[2]] + RotatedAxis())
  subP6 <- (pList4[[1]] + RotatedAxis()) + (pList4[[2]] + RotatedAxis())
  subP7 <- (pList5[[1]] + RotatedAxis()) + (pList5[[2]] + RotatedAxis())
  # Visualization key characteristics before and after filtering
  plotList <- list(subP1, subP2, subP3, subP4, subP5, subP6, subP7)
  do.call(gridExtra::grid.arrange, c(plotList, ncol=3))
  bafiltP <- do.call(gridExtra::arrangeGrob, c(plotList, ncol=2))
  bafiltP
  ggsave(filename = "before_and_after_basic_filtering.pdf", plot = bafiltP, device="pdf", width = 24, height = 18)
  
  stopCluster(cl)
  
  rm(list=c('sceList.old', 'sceList.merged'))
  # saveRDS(sceList, file = "sceList.Rds")
}

# Plot cell cycle details
# S.Score较高的为S期，G2M.Score较高的为G2M期，都比较低的为G1期
if(F){
  ccList <- list()
  for (i in seq(length(sceList))) {
    fname <- paste0("ccP",i)
    fname <- sceList[[i]]@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+ theme_minimal() + ggtitle(samples[i])
    ccList[[i]] <- fname
  }
  
  do.call(gridExtra::grid.arrange, c(ccList, ncol=2))
  ccP <- do.call(gridExtra::arrangeGrob, c(ccList, ncol=2))
  ccP
  ggsave(filename="cycle_details.pdf", plot = ccP, device = "pdf", width = 14, height = 10)
}

# Plot top 50 most expressed genes.
# May need to modify some codes if you have more than two samples.
if(F) {
  for (i in 1:length(sceList)) {
    C=sceList[[i]]@assays$SCT@counts
    dim(C)
    C=Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
    # If you don't have enough memory or cores to calculate, consider sampling here
    # C=C[,sample(1:ncol(C),1000)]
    most_expressed_C <- order(apply(C, 1, median), decreasing = T)[50:1]
    
    fname <- paste0("TOP50_most_expressed_gene_in_",samples[i],".pdf")
    pdf(file=fname, width = 16, height = 10)
    boxplot(as.matrix(Matrix::t(C[most_expressed_C, ])),
            cex = 0.1, las = 1, 
            xlab = "% total count per cell", 
            col = (scales::hue_pal())(50)[50:1], 
            horizontal = TRUE,
            main = strsplit(fname, ".",fixed=T)[[1]][1])
    dev.off()
    rm(C)
  }
}

if(F) {
  sceList <- lapply(X = sceList, FUN = function(x) {
    x <- RunPCA(x,npcs = 30, verbose = F) %>%
      RunUMAP(reduction = "pca", dims = 1:30, verbose = F) %>%
      FindNeighbors(dims = 1:30, verbose = FALSE) %>%
      FindClusters(verbose = FALSE)
  })
}
#################### Step 3. Find Doublets ####################
# https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/77
# Doubletfinder to remove doublets
if(T) {
  sweep.stats.list <- list()
  pk.vec <- vector()
  for (i in 1:length(sceList)) {
    seu_temp <- sceList[[i]]
    sweep.res.list <- paramSweep_v3(seu_temp, PCs = seu_temp@commands$RunUMAP.SCT.pca$dims, sct = TRUE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    sweep.stats.list[[i]] <- sweep.stats
    # You can then manually select your pK values for each Seurat object and store them in vector called pk.vec before 
    # iterating through the DoubletFinder run (assuming doublet formation rate of 7.5%)
    bcmvn <- find.pK(sweep.stats.list[[i]])
    pk.vec[i] <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  }
  
  df.list <- vector()
  # assign doublet formation rate according to 10X genomics estimations
  # https://uofuhealth.utah.edu/huntsman/shared-resources/gba/htg/single-cell/genomics-10x.php
  for (i in 1:length(sceList)) {
    seu_temp <- sceList[[i]]
    if (ncol(seu_temp) > 10000) {
      doublet_rate <- 0.076
    } else if (ncol(seu_temp) > 9000) {
      doublet_rate <- 0.069
    } else if (ncol(seu_temp) > 8000) {
      doublet_rate <- 0.061
    } else if (ncol(seu_temp) > 7000) {
      doublet_rate <- 0.054
    } else if (ncol(seu_temp) > 6000) {
      doublet_rate <- 0.046
    } else if (ncol(seu_temp) > 5000) {
      doublet_rate <- 0.039
    } else if (ncol(seu_temp) > 4000) {
      doublet_rate <- 0.031
    } else if (ncol(seu_temp) > 3000) {
      doublet_rate <- 0.023
    } else if (ncol(seu_temp) > 2000) {
      doublet_rate <- 0.016
    } else if (ncol(seu_temp) > 1000) {
      doublet_rate <- 0.008
    } else if (ncol(seu_temp) > 500) {
      doublet_rate <- 0.004
    } else {
      doublet_rate <- 0
    }
    
    nExp_poi <- doublet_rate*nrow(seu_temp@meta.data)
    df.list[i] <- paste0("DF.classifications_0.25_",pk.vec[i],"_",nExp_poi)
    
    seu_temp <- doubletFinder_v3(seu_temp, PCs = seu_temp@commands$RunUMAP.SCT.pca$dims, pN = 0.25,
                                 pK = pk.vec[i], nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
    
    lo = grep(df.list[i],names(seu_temp@meta.data))
    doublets <- seu_temp@meta.data[,c(lo-1, lo)]
    colnames(doublets) <- c("Doublet_score","Is_doublet")
    head(doublets)
    seu_temp <- AddMetaData(seu_temp,doublets)
    seu_temp[['QC']] <- ifelse(seu_temp@meta.data$Is_doublet == 'Doublet','Doublet','Pass')
    seu_temp[['QC']] <- ifelse(seu_temp@meta.data$nFeature_RNA < 200 & seu_temp@meta.data$QC == 'Pass','Low_nFeature',seu_temp@meta.data$QC)
    seu_temp[['QC']] <- ifelse(seu_temp@meta.data$nFeature_RNA < 200 & seu_temp@meta.data$QC != 'Pass' & seu_temp@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',seu_temp@meta.data$QC,sep = ','),seu_temp@meta.data$QC)
    seu_temp[['QC']] <- ifelse(seu_temp@meta.data$percent.mito > 20 & seu_temp@meta.data$QC == 'Pass','High_MT',seu_temp@meta.data$QC)
    seu_temp[['QC']] <- ifelse(seu_temp@meta.data$nFeature_RNA < 200 & seu_temp@meta.data$QC != 'Pass' & seu_temp@meta.data$QC != 'High_MT',paste('High_MT',seu_temp@meta.data$QC,sep = ','),seu_temp@meta.data$QC)
    table(seu_temp[['QC']])
    seu_temp <- subset(seu_temp, subset = QC == "Pass")
    sceList[[i]] <- seu_temp
    rm(seu_temp)
  }
  saveRDS(sceList, "sceList.no.doublets.Rds")
}

#################### Step 4. Data Integration and/or Batch correction using Seurat, harmony or SCTransformed based Integration ####################
# Choose integration method here
# 0. Direct merge multiple samples
# 1. Integration using Seurat
# 2. Batch correction using Harmony
# 3. Performing integration on datasets normalized with SCTransform
if(file.exists("sceList.no.doublets.Rds")) {
  #sceList <- readRDS("sceList.no.doublets.Rds")
}
# Choose from 0 to 3, default value is 1
intOpt <- 2

# ******** Option 0 (Default: off) ********* 
# Option 0: Direct merge multiple datasets
if(intOpt==0) {
  # Direct merge multiple datasets
  for (i in 1:length(sceList)) {
    DefaultAssay(sceList[[i]]) <- "RNA"
  }
  sce.all <- merge(sceList[[1]], y= sceList[-1] ,add.cell.ids = samples) 
  
  sce.all <- NormalizeData(sce.all) %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30) %>%
    RunTSNE(reduction = "pca", dims = 1:30) %>%
    FindNeighbors(reduction = "pca", k.param = 20, dims = 1:30) %>% 
    FindClusters(resolution = 0.8) %>% 
    identity()  
  
  umapP0 <- DimPlot(sce.all, reduction = "umap", split.by = 'orig.ident', label = T) + plot_annotation(title = "Umap plot of samples by directly merging multiple datasets")
  tsneP0 <- DimPlot(sce.all, reduction = "tsne", split.by = 'orig.ident', label = T) + plot_annotation(title = "Tsne plot of samples by directly merging multiple datasets")
  
  library(patchwork)
  p <- umapP0 + tsneP0
  ggsave(filename = "umap_tsne_plots_by_directly_merging_multiple_datasets.pdf", plot = p, device = "pdf", width = 24, height = 12)
  
  saveRDS(sce.all,file="sce.all.raw.Rds")
  
}

# ******** Option 1 (Default: on) ********* 
# Option 1: Integration using Seurat

if(intOpt==1) {
  # normalize and identify variable features for each dataset independently
  for (i in seq(length(sceList))) {
    DefaultAssay(sceList[[i]]) <- "RNA"
    sceList[[i]] <- NormalizeData(sceList[[i]]) %>% 
      FindVariableFeatures(selection.method = "vst", nfeatures = 3000)
  }
  
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = sceList)
  anchors <- FindIntegrationAnchors(object.list = sceList, anchor.features = features)
  sce.all <- IntegrateData(anchorset = anchors, new.assay.name = "CCA")
  
  DefaultAssay(sce.all) # "CCA"
  
  sce.all <- ScaleData(sce.all, verbose = F) %>%
    RunPCA(npcs = 30, verbose = F) %>%
    RunUMAP(reduction = "pca", dims = 1:30, verbose = F) %>%
    RunTSNE(reduction = "pca", dims = 1:30, verbose = F) %>%
    FindNeighbors(reduction = "pca", k.param = 20, dims = 1:30) %>% 
    FindClusters(resolution = 0.8) %>% 
    identity()
  
  # Plot UMAP and TSNE plot after integration using Seurat ingegration
  umapP1 <- DimPlot(sce.all, reduction = "umap", split.by = 'orig.ident', label = T) + plot_annotation(title = "Umap plot of samples after Seurat integration")
  tsneP1 <- DimPlot(sce.all, reduction = "tsne", split.by = 'orig.ident', label = T) + plot_annotation(title = "Tsne plot of samples after Seurat integration")
  
  library(patchwork)
  umapP1 + tsneP1
  
  umapP1
  ggsave(filename = "umap_plots_after_using_Seurat_integration.pdf", plot = umapP1, device = "pdf", width = 24, height = 12)
  tsneP1
  ggsave(filename = "tsne_plots_after_using_Seurat_integration.pdf", plot = tsneP1, device = "pdf", width = 24, height = 12)
  
  saveRDS(sce.all, "sce.all.seurat.Rds")
  
}

# ******** Option 2 (Default: off) ********* 
## Batch correction using Harmony
## batch correction can be used by using Seuart scRNA-seq intergration tools or Harmony
## Use either one should be fine.
# Ref: https://www.singlecellcourse.org/scrna-seq-dataset-integration.html

if(intOpt==2) {
  # Batch correction using harmony
  for (i in seq(length(sceList))) {
    DefaultAssay(sceList[[i]]) <- "RNA"
    sceList[[i]] <- NormalizeData(sceList[[i]]) %>% 
      FindVariableFeatures(selection.method = "vst", nfeatures = 3000)
  }
  sce.all <- merge(sceList[[1]], y= sceList[-1] ,add.cell.ids = samples) 
  DefaultAssay(sce.all) <- "RNA"
  
  sce.all <- NormalizeData(sce.all) %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30) %>%
    RunTSNE(reduction = "pca", dims = 1:30) %>%
    FindNeighbors(reduction = "pca", k.param = 20, dims = 1:30) %>% 
    FindClusters(resolution = 0.8) %>% 
    identity()
  
  umapP20 <- DimPlot(sce.all, reduction = "umap", split.by = 'orig.ident', label = T) + plot_annotation(title = "Umap plot of samples before harmony integration")
  tsneP20 <- DimPlot(sce.all, reduction = "tsne", split.by = 'orig.ident', label = T) + plot_annotation(title = "Tsne plot of samples before harmony integration")
  
  library(patchwork)
  p20 <- umapP20 + tsneP20
  ggsave(filename = "umap_tsne_plots_before_using_harmony.pdf", plot = p20, device = "pdf", width = 24, height = 12)
  
  sce.all <- sce.all %>% RunHarmony(group.by.vars = "orig.ident", assay.use="RNA", reduction = "pca", plot_convergence = TRUE)
  harmony_embeddings <- Embeddings(sce.all, 'harmony')
  harmony_embeddings[1:5, 1:5]
  
  sce.all <- sce.all %>% 
    RunUMAP(reduction = "harmony", dims = 1:30, verbose = F) %>% 
    RunTSNE(reduction = "harmony", dims = 1:30, verbose = F) %>% 
    FindNeighbors(reduction = "harmony", k.param = 20, dims = 1:30) %>% 
    FindClusters(resolution = 0.8) %>% 
    identity()
  
  # Plot UMAP and TSNE plot after integration using harmony
  umapP2 <- DimPlot(sce.all, reduction = "umap", split.by = 'orig.ident', label = T) + plot_annotation(title = "Umap plot of samples after harmony integration")
  tsneP2 <- DimPlot(sce.all, reduction = "tsne", split.by = 'orig.ident', label = T) + plot_annotation(title = "Tsne plot of samples after harmony integration")
  
  library(patchwork)
  p21 <- umapP2 + tsneP2
  ggsave(filename = "umap_tsne_plots_after_using_harmony.pdf", plot = p21, device = "pdf", width = 24, height = 12)
  
  saveRDS(sce.all,file="sce.all.harmony.Rds")
}

# ******** Option 3 (Default: off) ********* 
# Option 3: Performing integration on datasets normalized with SCTransform
# https://satijalab.org/seurat/articles/integration_introduction.html
# Integration on datasets using Seurat 

if(intOpt==3) {
  for (i in 1:length(sceList)) {
    DefaultAssay(sceList[[i]]) < "RNA"
  }
  sceList <- lapply(sceList, FUN = SCTransform)
  # Performing integration on datasets normalized with SCTransform
  features <- SelectIntegrationFeatures(object.list = sceList, nfeatures = 3000)
  sceList <- PrepSCTIntegration(object.list = sceList, anchor.features = features)
  
  anchors <- FindIntegrationAnchors(object.list = sceList, normalization.method = "SCT", anchor.features = features)
  sce.all <- IntegrateData(anchorset = anchors, new.assay.name = "CCA", normalization.method = "SCT")
  
  DefaultAssay(sce.all) <- "CCA"
  
  sce.all <- RunPCA(sce.all, npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30) %>%
    RunTSNE(reduction = "pca", dims = 1:30) %>%
    FindNeighbors(reduction = "pca", k.param = 20, dims = 1:30) %>% 
    FindClusters(resolution = 0.8) %>% 
    identity()
  
  sce.all <- RunUMAP(sce.all, reduction = "pca", dims = 1:30)
  sce.all <- RunTSNE(sce.all, reduction = "pca", dims = 1:30)
  
  saveRDS(sce.all, file = "sce.all.Rds")
  # Plot UMAP and TSNE plot after integration using Seurat
  umapP3 <- DimPlot(sce.all, reduction = "umap", split.by = 'orig.ident', label = T) + plot_annotation(title = "Umap plot of samples after SCTransform integration")
  tsneP3 <- DimPlot(sce.all, reduction = "tsne", split.by = 'orig.ident', label = T) + plot_annotation(title = "Tsne plot of samples after SCTransform integration")
  
  library(patchwork)
  p <- umapP3 + tsneP3
  ggsave(filename = "umap_tsne_plots_after_SCTransform_integration.pdf", plot = p, device = "pdf", width = 24, height = 12)
  
  saveRDS(sce.all, "sce.all.sct.Rds")
}


#################### Step 5. Find best dimmensions ####################
## Find best dimmensions.
## This is jsut for reference.

if(file.exists("sce.all.seurat.Rds")) {
  sce.all <- readRDS("sce.all.seurat.Rds")
}

if(T) {
  # set random seed to make duplicable results
  set.seed(010101)
  stdev <- sce.all@reductions$pca@stdev
  var <- stdev^2
  EndVar = 0
  for(i in 1:length(var)){
    total <- sum(var)
    numerator <- sum(var[1:i])
    expvar <- numerator/total
    if(EndVar == 0){
      if(expvar > 0.90){
        EndVar <- EndVar + 1
        PCNum <- i
      }
    }
  }
  
  ## Confirm PC's determined explain > 90% of variance
  sum(var[1:PCNum])/sum(var)
  print(PCNum)
  
  visP <- VizDimLoadings(sce.all, dims = 1:12, reduction = "pca")
  ggsave(filename = "vizdimloading_of_sce.all_dim1_2.pdf", plot = visP, device = "pdf", width = 12, height = 8)
  
  sce.all <- RunUMAP(sce.all, reduction = "harmony", dims = 1:PCNum)
  sce.all <- RunTSNE(sce.all, reduction = "harmony", dims = 1:PCNum)
}

#################### Step 7. Find best resolution ####################
# The result is only for reference
# Set different resolutions, choose the best one for downstream analysis.
if(T) {
  library(clustree)
  sce.all <- FindNeighbors(sce.all, dims = 1:PCNum, k.param = 20, prune.SNN = 1/15)
  sce.all <- FindClusters(object = sce.all, resolution = c(seq(0.1,2.0,0.1)))
  
  library(cowplot)
  p1_dim=plot_grid(DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.1", label = T, raster=F) + ggtitle("louvain_0.1"),
                   DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.3", label = T, raster=F) + ggtitle("louvain_0.3"),
                   DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.5", label = T, raster=F) + ggtitle("louvain_0.5"),
                   DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.0.8", label = T, raster=F) + ggtitle("louvain_0.8"),
                   DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.1", label = T, raster=F) + ggtitle("louvain_1"),
                   DimPlot(sce.all, reduction = "umap", group.by = "RNA_snn_res.1.2", label = T, raster=F) + ggtitle("louvain_1.2"),
                   ncol = 3)
  ggsave(filename="Dimplot_6diff_resolution.pdf", plot=p1_dim, device = "pdf", width = 28, height = 20)
  
  ctreeP <- clustree(sce.all@meta.data, prefix = "RNA_snn_res.")
  ggsave(filename = "clustree_diff_resolutions.pdf", plot =ctreeP, device = "pdf", width = 10, height = 14)
}

# You need to check UMAP or TSNE plot to determine the best resolution
# Set the best resolution here, default values is 0.3
bestRes <- 0.8
if(bestRes > 0) {
  sce.all <- FindNeighbors(sce.all, dims = 1:PCNum) %>%
    FindClusters(resolution = bestRes)
  
  gpb <- paste0("RNA_snn_res.",bestRes)
  bestDimUMAPPlot <- DimPlot(sce.all, reduction = "umap", group.by = gpb, split.by = 'orig.ident', label = T, raster=F) + ggtitle(gpb)
  ggsave(filename = "best_dim_umap_plot_sce.all.pdf", plot = bestDimUMAPPlot, device = "pdf", width = 14, height = 8)
  
  bestDimTSNEPlot <- DimPlot(sce.all, reduction = "tsne", group.by = gpb, split.by = 'orig.ident', label = T, raster=F) + ggtitle(gpb)
  ggsave(filename = "best_dim_tsne_plot_sce.all.pdf", plot = bestDimTSNEPlot, device = "pdf", width = 14, height = 8)
  
  saveRDS(sce.all, "sce.all_best_resolution.Rds")
}

# Visualize the distribution of cells of different datasets per cluster, alongside cluster sizes:
# https://www.singlecellcourse.org/scrna-seq-dataset-integration.html
if(T) {
  picP <- plot_integrated_clusters(sce.all)
  ggsave(file = "distribtion_of_cells.pdf",plot = picP, device = "pdf", width = 12, height = 8)
}

source("../step2-CheckCellMarkers.R")