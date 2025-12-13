########################### Monocle2 for Cd4 T cells ######################
setwd("/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype/newTcells/")

dir.create("Cd4")
setwd("Cd4")
sce.t <- readRDS("/home/LiuLab/zzdx/data/singlecell/bgi/wangpengju/xuanyujing/results/4-modifiedcelltype/sce.T_w_correction.Rds")

Idents(sce.t) <- sce.t$newsubcelltype
sce.t.cd4 <- subset(sce.t, idents = c('Activated Th2', 'Cd4 T early activated', 'Cd4 Tcm', 'Cd4 Tem', 
                                        'Proliferating Treg','Tfh','Th1','Th2','Tn','Treg'))

saveRDS(sce.t.cd4, "sce.t.Cd4.Rds")

sce.t.cd8 <- subset(sce.t, idents = c('Cd8 Teff', 'Cd8 Tex', 'Cd8 Tm', 'Cd8 Tpex', 'IFN-responsive Cd8 Teff','Ly6c+ Cd8 Teff',
                                      'Proliferating Cd8 T','Proliferating Cd8 Teff','Tn'))

saveRDS(sce.t.cd8, "sce.t.Cd8.Rds")


## ----------------------
## 0. 准备：从 Seurat 到 Monocle2
## ----------------------
library(Seurat)
library(monocle)

# 建议固定随机种子，方便复现轨迹形状
set.seed(1234)

# 用更显式的方式指定 celltype
sce.t.cd4$celltype <- sce.t.cd4$newsubcelltype
Idents(sce.t.cd4) <- "celltype"
table(Idents(sce.t.cd4))

seurat <- sce.t.cd4

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

# 只取前 2000 个最显著的基因（按需要调整）
deg_filtered <- deg_filtered[order(deg_filtered$qval), ]
ordering_genes <- rownames(head(deg_filtered, 2000))

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
saveRDS(cds, "cds.cd4_monocle2.rds")



## ----------------------
## 4. 细胞排序：指定 root_state（如已知起点群）
## ----------------------

# 先看一下 celltype/cluster 在轨迹上的分布
plot_cell_trajectory(cds, color_by = "celltype")
plot_cell_trajectory(cds, color_by = "orig.ident")
plot_cell_trajectory(cds, color_by = "celltype") + facet_wrap("~celltype", nrow = 2)
plot_cell_trajectory(cds, color_by = "group") + facet_wrap("~group", nrow = 2)

p <- plot_complex_cell_trajectory(cds = cds, x=1,y=2, color_by = "celltype") +
  theme(legend.title = element_blank()) 
p
cds <- orderCells(cds)


p <- plot_complex_cell_trajectory(cds = cds.cd8, x=1,y=2, color_by = "celltype") +
  theme(legend.title = element_blank()) 
p

plot_cell_trajectory(cds.cd8, color_by = "celltype") + facet_wrap("~celltype", nrow = 2)
# 如果你已经知道某个状态是“起点”，可以手动设 root_state
# 例：假设我们想让包含 naive CD4 T 的 branch 当起点
# 先看看每个 state 中 celltype 的组成：
state_tab <- table(pData(cds)$State, pData(cds)$celltype)

# 确认一下 Tn 这一列是否存在
colnames(state_tab)

# 找到 Tn 数量最多的那个 State
state_for_Tn <- rownames(state_tab)[which.max(state_tab[, "Tn"])]
state_for_Tn

# 假设我们选 state 3 为起点（根据上面的表判断）
root_state_to_use <- 8

cds <- orderCells(
  cds,
  root_state = root_state_to_use
)

# 如果暂时不确定起点，也可以先不填 root_state：
# cds <- orderCells(cds)

## ----------------------
## 5. 可视化：多种颜色 / 分面
## ----------------------

# 1) 按 orig.ident 着色
p1 <- plot_cell_trajectory(cds, color_by = "orig.ident")

# 2) 按 celltype 着色
p2 <- plot_cell_trajectory(cds, color_by = "celltype")

# 3) 分面看每个 celltype 在轨迹上的位置
p3 <- plot_cell_trajectory(cds, color_by = "celltype") +
  facet_wrap(~ celltype, nrow = 2)

# 4) 按伪时间上色
p4 <- plot_cell_trajectory(cds, color_by = "Pseudotime")

# 5) 按 seurat_clusters 着色（前提是 meta.data 里有这个列）
p5 <- plot_cell_trajectory(cds, color_by = "seurat_clusters")

# 6) 调整点大小
p6 <- plot_cell_trajectory(cds, color_by = "celltype", cell_size = 0.8)

# 按需打印或保存
print(p1); print(p2); print(p3); print(p4); print(p5); print(p6)

# 自动找到“Tn 数量最多的 State”
state_tab <- table(pData(cds)$State, pData(cds)$celltype)

# 确认一下 Tn 这一列是否存在
colnames(state_tab)

# 找到 Tn 数量最多的那个 State
state_for_Tn <- rownames(state_tab)[which.max(state_tab[, "Tn"])]
state_for_Tn

cds <- orderCells(cds, root_state = as.numeric(state_for_Tn))


plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 2)

plot_cell_trajectory(cds,cell_size = 1)

plot_cell_trajectory(cds,cell_size = 1, color_by = "Pseudotime")
plot_cell_trajectory(cds,cell_size = 1, color_by = "seurat_clusters")

sce <- SCTransform(scRNA,verbose = F)  
sce <- RunPCA(sce,verbose = F)  
sce <- RunUMAP(sce, dims = 1:20)
DimPlot(sce,group.by = 'celltype')

sce$State = pData(cds)$State
sce$Pseudotime = pData(cds)$Pseudotime
DimPlot(sce,group.by = 'State', label = T) + 
  DimPlot(sce,group.by = 'celltype', label = T) +
  DimPlot(sce,group.by = 'orig.ident', label = T) 

