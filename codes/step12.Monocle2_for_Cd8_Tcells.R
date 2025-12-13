########################### Monocle2 for Cd4 T cells ######################
setwd("/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype/newTcells/Cd8")

#sce.t <- readRDS("/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype/sce.T_w_correction.Rds")

#Idents(sce.t) <- sce.t$newsubcelltype
#sce.t.cd4 <- subset(sce.t, idents = c('Activated Th2', 'Cd4 T early activated', 'Cd4 Tcm', 'Cd4 Tem', 
#                                        'Proliferating Treg','Tfh','Th1','Th2','Tn','Treg'))

#saveRDS(sce.t.cd4, "sce.t.Cd4.Rds")

#sce.t.cd8 <- subset(sce.t, idents = c('Cd8 Teff', 'Cd8 Tex', 'Cd8 Tm', 'Cd8 Tpex', 'IFN-responsive Cd8 Teff','Ly6c+ Cd8 Teff',
#                                      'Proliferating Cd8 T','Proliferating Cd8 Teff','Tn'))

#saveRDS(sce.t.cd8, "sce.t.Cd8.Rds")
sce.t.cd8 <- readRDS("sce.t.Cd8.Rds")
## ----------------------
## 0. 准备：从 Seurat 到 Monocle2
## ----------------------
library(Seurat)
library(monocle)

# 建议固定随机种子，方便复现轨迹形状
set.seed(1234)

# 用更显式的方式指定 celltype
sce.t.cd8$celltype <- sce.t.cd8$newsubcelltype
Idents(sce.t.cd8) <- "celltype"
table(Idents(sce.t.cd8))

seurat <- sce.t.cd8

# 确保使用 RNA assay（如果你对象里有多个 assay）
DefaultAssay(seurat) <- "RNA"

## ----------------------
## 1. 构建 Monocle2 的 cds 对象
## ----------------------

# 建议使用 counts（稀疏矩阵），适合 negbinomial
expr_matrix <- GetAssayData(seurat, assay = "RNA", slot = "counts")

# phenotype data
sample_sheet <- seurat@meta.data

# feature data：建议加上 gene_id，方便后续注释
gene_annotation <- data.frame(
  gene_short_name = rownames(seurat),
  row.names = rownames(seurat)
)

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

cds <- newCellDataSet(
  expr_matrix,
  phenoData  = pd,
  featureData = fd,
  expressionFamily = VGAM::negbinomial.size()
)

# 估计 size factor / dispersion
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

## ----------------------
## 2. 选择 ordering genes（关键！）
##    用差异分析 + 简单过滤，避免垃圾基因主导轨迹
## ----------------------

# 差异分析：考虑 celltype，同时把 orig.ident 当作协变量
diff_test_res <- differentialGeneTest(
  cds,
  fullModelFormulaStr    = "~ celltype + orig.ident",
  reducedModelFormulaStr = "~ orig.ident",
  relative_expr = TRUE,
  cores = 8
)

# 按 qval 过滤，再限制一下基因数，避免过多
deg <- subset(diff_test_res, qval < 0.01)

# 可选：去掉线粒体、核糖体等技术噪音基因（按需开启）
rm_mt <- grepl("^mt-", rownames(deg), ignore.case = TRUE)
rm_ribo <- grepl("^Rpl|^Rps", rownames(deg), ignore.case = TRUE)
deg_filtered <- deg[!(rm_mt | rm_ribo), ]

# 只取前 1000 个最显著的基因（按需要调整）
deg_filtered <- deg_filtered[order(deg_filtered$qval), ]
ordering_genes <- rownames(head(deg_filtered, 1000))

length(ordering_genes)
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)

## ----------------------
## 3. 降维：DDRTree 参数优化
## ----------------------

cds <- reduceDimension(
  cds,
  max_components = 2,
  num_dim = 20,               # 稍微多一点 PCA 维度，轨迹更稳定
  reduction_method = "DDRTree",
  residualModelFormulaStr = "~ orig.ident",  # 去掉样本批次影响
  norm_method = "log",        # 强烈建议，避免 count 尺度过大
  verbose = TRUE
)

# 保存一份中间结果
saveRDS(cds, "cds.cd3_monocle2.rds")
