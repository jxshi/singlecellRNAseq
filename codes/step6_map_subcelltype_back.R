setwd("/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/3-celltype/Tcells")

sce.epcells <- readRDS("../Epithelialcells/Epithelialcells.Rds")
sce.epcells$subcelltype <- "Epithelial cells"
sce.epcells$maincelltype <- "Epithelial cells"

sce.fibro <- readRDS("../Fibroblasts/Fibroblasts.Rds")
sce.fibro$subcelltype <- "Fibroblasts"
sce.fibro$maincelltype <- "Fibroblasts"

sce.gran <- readRDS("../Neutrophils/Neutrophils.Rds")
sce.gran$subcelltype <- "Granulocytes"
sce.gran$maincelltype <- "Granulocytes"

setwd("/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype/")
sce.all <- readRDS("../3-celltype/sce.all_w_cell_type_anno.Rds")

# Add the subcelltype to the original Seurat object metadata
sce.all@meta.data$bsubcelltype <- sce.b@meta.data$newsubcelltype[match(rownames(sce.all@meta.data), rownames(sce.b@meta.data))]
sce.all@meta.data$tsubcelltype <- sce.t@meta.data$newsubcelltype[match(rownames(sce.all@meta.data), rownames(sce.t@meta.data))]
sce.all@meta.data$dcsubcelltype <- sce.dc@meta.data$newsubcelltype[match(rownames(sce.all@meta.data), rownames(sce.dc@meta.data))]
sce.all@meta.data$monosubcelltype <- sce.monomacro@meta.data$newsubcelltype[match(rownames(sce.all@meta.data), rownames(sce.monomacro@meta.data))]
sce.all@meta.data$gransubcelltype <- sce.gran@meta.data$subcelltype[match(rownames(sce.all@meta.data), rownames(sce.gran@meta.data))]
sce.all@meta.data$fibrosubcelltype <- sce.fibro@meta.data$subcelltype[match(rownames(sce.all@meta.data), rownames(sce.fibro@meta.data))]
sce.all@meta.data$episubcelltype <- sce.epcells@meta.data$subcelltype[match(rownames(sce.all@meta.data), rownames(sce.epcells@meta.data))]
sce.all@meta.data$mastsubcelltype <- sce.mast@meta.data$subcelltype[match(rownames(sce.all@meta.data), rownames(sce.mast@meta.data))]


table(sce.all$bsubcelltype, sce.all$group)
table(sce.all$tsubcelltype, sce.all$group)
table(sce.all$dcsubcelltype, sce.all$group)
table(sce.all$monosubcelltype, sce.all$group)
table(sce.all$gransubcelltype, sce.all$group)
table(sce.all$fibrosubcelltype, sce.all$group)
table(sce.all$episubcelltype, sce.all$group)
table(sce.all$mastsubcelltype, sce.all$group)


table(sce.all@meta.data$celltype)

library(dplyr)
sce.all$newsubcelltype <- ""
subclty <- sce.all@meta.data

subclty <- subclty %>%
  mutate(subcelltype = case_when(
    !is.na(bsubcelltype) & bsubcelltype != "" ~ bsubcelltype,
    !is.na(tsubcelltype) & tsubcelltype != "" ~ tsubcelltype,
    !is.na(dcsubcelltype) & dcsubcelltype != "" ~ dcsubcelltype,
    !is.na(monosubcelltype) & monosubcelltype != "" ~ monosubcelltype,
    !is.na(gransubcelltype) & gransubcelltype != "" ~ gransubcelltype,
    !is.na(fibrosubcelltype) & fibrosubcelltype != "" ~ fibrosubcelltype,
    !is.na(episubcelltype) & episubcelltype != "" ~ episubcelltype,
    !is.na(mastsubcelltype) & mastsubcelltype != "" ~ mastsubcelltype,
    !is.na(subcelltype) & subcelltype != "" ~ subcelltype,
    TRUE ~ NA_character_  # Assign NA to subcelltype if none of the conditions are met
  ))
table(subclty$subcelltype)
dim(subclty)

length(sce.all$subcelltype)

sce.all$subcelltype <- subclty$subcelltype

table(sce.all$celltype, sce.all$subcelltype)


######################################## Add maincelltype back to original sce.all object ##################################
# Add the modified main celltype to the original Seurat object metadata
sce.all@meta.data$bmaincelltype <- sce.b@meta.data$maincelltype[match(rownames(sce.all@meta.data), rownames(sce.b@meta.data))]
sce.all@meta.data$tmaincelltype <- sce.t@meta.data$maincelltype[match(rownames(sce.all@meta.data), rownames(sce.t@meta.data))]
sce.all@meta.data$dcmaincelltype <- sce.dc@meta.data$maincelltype[match(rownames(sce.all@meta.data), rownames(sce.dc@meta.data))]
sce.all@meta.data$monomaincelltype <- sce.monomacro@meta.data$maincelltype[match(rownames(sce.all@meta.data), rownames(sce.monomacro@meta.data))]
sce.all@meta.data$granmaincelltype <- sce.gran@meta.data$maincelltype[match(rownames(sce.all@meta.data), rownames(sce.gran@meta.data))]
sce.all@meta.data$fibromaincelltype <- sce.fibro@meta.data$maincelltype[match(rownames(sce.all@meta.data), rownames(sce.fibro@meta.data))]
sce.all@meta.data$epimaincelltype <- sce.epcells@meta.data$maincelltype[match(rownames(sce.all@meta.data), rownames(sce.epcells@meta.data))]
sce.all@meta.data$mastmaincelltype <- sce.mast@meta.data$maincelltype[match(rownames(sce.all@meta.data), rownames(sce.mast@meta.data))]


table(sce.all$bmaincelltype, sce.all$group)
table(sce.all$tmaincelltype, sce.all$group)
table(sce.all$dcmaincelltype, sce.all$group)
table(sce.all$monomaincelltype, sce.all$group)
table(sce.all$granmaincelltype, sce.all$group)
table(sce.all$fibromaincelltype, sce.all$group)
table(sce.all$epimaincelltype, sce.all$group)
table(sce.all$mastmaincelltype, sce.all$group)


table(sce.all@meta.data$celltype)

library(dplyr)
sce.all$maincelltype <- ""
subclty <- sce.all@meta.data

subclty <- subclty %>%
  mutate(maincelltype = case_when(
    !is.na(bmaincelltype) & bmaincelltype != "" ~ bmaincelltype,
    !is.na(tmaincelltype) & tmaincelltype != "" ~ tmaincelltype,
    !is.na(dcmaincelltype) & dcmaincelltype != "" ~ dcmaincelltype,
    !is.na(monomaincelltype) & monomaincelltype != "" ~ monomaincelltype,
    !is.na(granmaincelltype) & granmaincelltype != "" ~ granmaincelltype,
    !is.na(fibromaincelltype) & fibromaincelltype != "" ~ fibromaincelltype,
    !is.na(epimaincelltype) & epimaincelltype != "" ~ epimaincelltype,
    !is.na(mastmaincelltype) & mastmaincelltype != "" ~ mastmaincelltype,
    !is.na(maincelltype) & maincelltype != "" ~ maincelltype,
    TRUE ~ NA_character_  # Assign NA to maincelltype if none of the conditions are met
  ))
table(subclty$maincelltype)
dim(subclty)

length(sce.all$maincelltype)

sce.all$maincelltype <- subclty$maincelltype

table(sce.all$celltype, sce.all$maincelltype)

DimPlot(sce.all, raster = F, label = T)
Idents(sce.all) <- sce.all$maincelltype
DimPlot(sce.all, raster = F, label = T)

saveRDS(sce.all, "sce.all_w_corrected_main_celltype.Rds")


sce.new <- subset(sce.all, subset = maincelltype != "Unknown")
DimPlot(sce.new, raster = F, label = T)

saveRDS(sce.new, "sce.new_w_corrected_main_celltype.Rds")
