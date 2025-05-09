# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("STdeconvolve")
remotes::install_github("JEFworks-Lab/STdeconvolve")
# install.packages("devtools")
devtools::install_github("ati-lz/ISCHIA")
# Loading required packages
library(ISCHIA)
library(robustbase)
library(data.table)
library(ggplot2)
library(Seurat)
library(dplyr)
library(stringr)
#### 测试数据 ####
# 提取去卷积的矩阵，可以是RCTD和Cell2location
load("IBD_visium_SeuratObj_small.RData")
deconv.mat <- as.matrix(IBD.visium.P4@meta.data[,9:28])
range(deconv.mat)
#0.0000000 0.7513638,说明是0-1标准化后的比例
colnames(deconv.mat) <- sapply(colnames(deconv.mat), function(x) unlist(strsplit(x, split = "_"))[2])

Composition.cluster.k(deconv.mat, 20)

# 使用 Composition.cluster 函数对解卷积（deconvoluted）的空间点进行聚类
IBD.visium.P4 <- Composition.cluster(IBD.visium.P4, deconv.mat, 8)
# 聚类结果不太理想，接下来将通过 "CompositionCluster_CC" 对象设置 Idents（标识）

# 检查每个聚类的点数分布
table(IBD.visium.P4$CompositionCluster_CC)
# 输出聚类结果：每个聚类（CC1 到 CC8）包含的点的数量
#> 
#> CC1 CC2 CC3 CC4 CC5 CC6 CC7 CC8 
#> 302 296 274 382 108 238 426 159

# 绘制 SpatialDimPlot，按聚类结果（CompositionCluster_CC）分组
# 并为每个聚类指定对应的颜色
SpatialDimPlot(IBD.visium.P4, group.by = c("CompositionCluster_CC")) +
  scale_fill_manual(values = c("cyan", "orange", "purple", "green", "yellow", "blue", "red", "black"))
#> Scale for fill is already present.
#> Adding another scale for fill, which will replace the existing scale.
# 检查每个 Composition Cluster (CC) 的主要细胞类型
Composition_cluster_enrichedCelltypes(IBD.visium.P4, "CC4", deconv.mat) 
# 检查聚类 CC4 的主要富集细胞类型
Composition_cluster_enrichedCelltypes(IBD.visium.P4, "CC7", deconv.mat) 
# 检查聚类 CC7 的主要富集细胞类型
# 探索基于解卷积细胞类型组成的 UMAP 表示
IBD.visium.P4.umap <- Composition_cluster_umap(IBD.visium.P4, deconv.mat) 
# 使用解卷积后的细胞类型组成生成空间点的 UMAP 表示
# 提示：如果数据集较大（如 2185 个像素，20 种细胞类型），此过程可能需要较长时间
# 可视化 UMAP 结果
IBD.visium.P4.umap$umap.cluster.gg 
# 绘制按 Composition Cluster (CC) 显示的 UMAP 图，颜色表示所属的 CC（如 CC1, CC2 等）
IBD.visium.P4.umap$umap.deconv.gg 
# 绘制按解卷积细胞类型组成显示的 UMAP 图
# 每个点用饼图表示，饼图各部分表示该空间点中不同细胞类型的比例
# 对 CC4 聚类进行细胞类型共现分析
CC4.celltype.cooccur <- spatial.celltype.cooccurence(
  spatial.object = IBD.visium.P4,        # 输入的空间数据对象
  deconv.prob.mat = deconv.mat,         # 解卷积的概率矩阵
  COI = "CC4",                          # 感兴趣的聚类 (Composition of Interest, COI)
  prob.th = 0.05,                       # 概率阈值，控制共现计算的显著性
  Condition = unique(IBD.visium.P4$orig.ident) # 条件变量（如样本来源）
)

# 打印共现分析摘要，并绘制 CC4 的细胞类型共现图
plot.celltype.cooccurence(CC4.celltype.cooccur)

# 对 CC7 聚类进行细胞类型共现分析
CC7.celltype.cooccur <- spatial.celltype.cooccurence(
  spatial.object = IBD.visium.P4,        # 输入的空间数据对象
  deconv.prob.mat = deconv.mat,         # 解卷积的概率矩阵
  COI = "CC7",                          # 感兴趣的聚类 (Composition of Interest, COI)
  prob.th = 0.05,                       # 概率阈值
  Condition = unique(IBD.visium.P4$orig.ident) # 条件变量（如样本来源）
)

# 打印共现分析摘要，并绘制 CC7 的细胞类型共现图
plot.celltype.cooccurence(CC7.celltype.cooccur)

# 由于配体-受体共现分析需要较长时间运行，目前将代码注释掉，仅保留逻辑

# 读取配体-受体网络数据
lr_network <- readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
all.LR.network <- cbind(
  lr_network[, c("from", "to")],                   # 配体和受体对
  LR_Pairs = paste(lr_network$from, lr_network$to, sep = "_") # 将配体和受体组合成字符串
)

# 筛选配体和受体基因存在于数据中的配体-受体对
all.LR.network.exp <- all.LR.network[
  which(all.LR.network$from %in% rownames(IBD.visium.P4) & 
          all.LR.network$to %in% rownames(IBD.visium.P4)),
]

# 为了减少计算时间，随机抽取 500 个配体-受体对作为样本
all.LR.network.exp <- sample_n(all.LR.network.exp, 500)

# 获取所有相关基因，并找到它们在数据中的交集
all.LR.genes <- unique(c(all.LR.network.exp$from, all.LR.network.exp$to))
all.LR.genes.comm <- intersect(all.LR.genes, rownames(IBD.visium.P4))

# 生成所有可能的配体-受体基因组合
LR.pairs <- all.LR.network.exp$LR_Pairs
LR.pairs.AllCombos <- combn(all.LR.genes.comm, 2, paste0, collapse = "_")

# 以下是具体的配体-受体共现分析代码（已注释掉以节约时间）

# CC4.Enriched.LRs <- Enriched.LRs(
#   IBD.visium.P4, c("CC4"), unique(IBD.visium.P4$orig.ident), 
#   all.LR.genes.comm, LR.pairs, 1, 0.2
# )

# CC7.Enriched.LRs <- Enriched.LRs(
#   IBD.visium.P4, c("CC7"), unique(IBD.visium.P4$orig.ident), 
#   all.LR.genes.comm, LR.pairs, 1, 0.2
# )

# 比较 CC4 和 CC7 中差异化显著的配体-受体对
# CC1vsCC4.Enriched.LRs.Specific <- Diff.cooc.LRs(
#   CC4.Enriched.LRs, CC7.Enriched.LRs, 0.05, 0.1
# )

# 绘制配体-受体共现结果的弦图
# ChordPlot.Enriched.LRs(CC4.Enriched.LRs$COI.enrcihed.LRs[1:20,])

# 绘制差异配体-受体对的桑基图
# SankeyPlot.Diff.LRs(CC4.Enriched.LRs$COI.enrcihed.LRs, CC7.Enriched.LRs$COI.enrcihed.LRs)

#### 实际运行 ####
# Loading required packages
library(ISCHIA)
library(robustbase)
library(data.table)
library(ggplot2)
library(Seurat)
library(dplyr)
library(factoextra)
library(cluster)
library(gridExtra)
library(tidyverse)
#install.packages("qpdf", type = "binary")
library(qpdf)
# 设置选项以避免科学计数法
options(scipen = 999)
# Set random seed for reproducibility
set.seed(123)

# Load data
load("G:/E-disk/ZJH_KRAS/Workdir/7.RCTD/ZJH_KRAS_RCTD_results/all_myRCTD_results.rda")
CRC_std=as.data.frame(rctd_weights)
# 检查标准化后的范围
range(CRC_std)
# Elbow Method
k.values <- 1:20
wss_values <- sapply(k.values, function(k) kmeans(CRC_std, k, nstart = 10)$tot.withinss)

pdf("1_elbow_plot.pdf")
plot(k.values, wss_values, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K", ylab = "Total within-cluster sum of squares",
     main = "Elbow Method for Optimal K")
dev.off()
elbow_k <- which.max(diff(wss_values))
# Gap Statistic
#Gap Statistic 方法的最佳 ????是第一个使 Gap 值显著增大并达到最大值的 
#对随机数据运行 k-均值聚类，计算离散度。重复 10 次
gap_stat <- function(k) {
  km.res <- kmeans(CRC_std, k, nstart = 10)
  if (k == 1) return(NA)
  obs_disp <- sum(km.res$withinss)
  reference_disp <- mean(replicate(10, {
    km.null <- kmeans(matrix(rnorm(nrow(CRC_std) * ncol(CRC_std)), 
                             ncol = ncol(CRC_std)), k, nstart = 10)
    sum(km.null$withinss)
  }))
  log(reference_disp) - log(obs_disp)
}
#下一步运行很慢
gap_stat_values <- sapply(k.values, gap_stat)

pdf("2_gap_statistic_plot.pdf")
plot(k.values, gap_stat_values, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of Clusters (K)", ylab = "Gap Statistic",
     main = "Gap Statistic: Determining Optimal K")
dev.off()

# Calinski-Harabasz Index 函数
calinski_harabasz_index <- function(data, labels) {
  num_clusters <- length(unique(labels))  # 簇的数量
  num_points <- nrow(data)                # 数据点总数
  
  # 计算全局中心点
  global_centroid <- colMeans(data)
  
  # 按簇分组，计算每个簇的中心点
  centroids <- sapply(unique(labels), function(cluster) {
    colMeans(data[labels == cluster, , drop = FALSE])  # 计算簇的中心点
  })
  centroids <- t(centroids)  # 转置矩阵
  
  # 计算簇间离散度
  between_disp <- sum(sapply(unique(labels), function(cluster) {
    cluster_size <- sum(labels == cluster)  # 簇中的数据点数量
    sum((centroids[cluster, ] - global_centroid) ^ 2) * cluster_size
  }))
  
  # 计算簇内离散度
  within_disp <- sum(sapply(unique(labels), function(cluster) {
    cluster_points <- data[labels == cluster, , drop = FALSE]
    sum(rowSums((cluster_points - centroids[cluster, ]) ^ 2))
  }))
  
  # Calinski-Harabasz 指数公式
  (between_disp / (num_clusters - 1)) / (within_disp / (num_points - num_clusters))
}

# 计算不同 k 值的 Calinski-Harabasz 指数
ch_values <- sapply(k.values, function(k) {
  km.res <- kmeans(CRC_std, k, nstart = 10)  # 对数据进行 K-means 聚类
  calinski_harabasz_index(CRC_std, km.res$cluster)  # 计算 Calinski-Harabasz 指数
})

# 绘制 Calinski-Harabasz 指数的图表
pdf("3_calinski_harabasz_plot.pdf")
plot(k.values, ch_values, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters (K)", ylab = "Calinski-Harabasz Index",
     main = "Calinski-Harabasz Index: Determining Optimal K")
dev.off()
# 找到 Calinski-Harabasz 指数的局部最大值位置


# ISCHIA Analysis
pdf("4_composition_cluster_k_plot.pdf")
Composition.cluster.k(CRC_std, 20)
dev.off()

ZJH_KRAS_combined_seurat <- readRDS("G:/E-disk/ZJH_KRAS/Workdir/6.Cluster/output/ZJH_KRAS_combined_seurat.rds")
samples=unique(ZJH_KRAS_combined_seurat@meta.data[["sample"]])
# 获取列名
rownames <- rownames(CRC_std)
# 循环修改列名
for (i in seq_along(samples)) {
  # 构造匹配的正则表达式
  pattern <- paste0(samples[i],"_")
  
  # 匹配列名并修改：添加前缀，删除 `_i`
  rownames <- ifelse(
    grepl(pattern, rownames),
    paste0(str_replace(rownames, pattern, ""),"_",i),
    rownames
  )
}

# 更新 CRC_sp 的列名
rownames(CRC_std) <- rownames
CRC_sp=ZJH_KRAS_combined_seurat[,rownames]

CRC_sp <- Composition.cluster(CRC_sp, CRC_std, 5)
CRC_sp$cc_5 <- CRC_sp$CompositionCluster_CC

# Spatial Dimension Plot
image_names <- samples

paletteMartin <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', 
                   '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff')

all_ccs <- unique(CRC_sp$CompositionCluster_CC)
color_mapping <- setNames(paletteMartin[1:length(all_ccs)], all_ccs)

pdf("5_spatial_plots_K5.pdf", width = 10, height = 7)
for (image_name in image_names) {
  plot <- SpatialDimPlot(CRC_sp, group.by = "CompositionCluster_CC", images = image_name) +
    scale_fill_manual(values = color_mapping) +
    theme_minimal() +
    ggtitle(image_name)
  print(plot)
}
dev.off()

# Enriched Cell Types
save_cc_plot <- function(cc) {
  plot <- Composition_cluster_enrichedCelltypes(CRC_sp, cc, as.matrix(CRC_std))
  pdf_name <- paste0(cc, ".pdf")
  pdf(file = pdf_name)
  print(plot)
  dev.off()
}

ccs <- paste0("CC", 1:5)
for (cc in ccs) {
  save_cc_plot(cc)
}

pdf_files <- paste0("CC", 1:5, ".pdf")
pdf_combine(pdf_files, output = "6_enrichedCelltypes_CC_5.pdf")

# UMAP
CRC_sp.umap <- Composition_cluster_umap(CRC_sp, CRC_std)
pdf("7_umap_pie_chart_CC_6.pdf")
print(CRC_sp.umap$umap.deconv.gg)
dev.off()

# Add UMAP to Seurat object
emb.umap <- CRC_sp.umap$umap.table
emb.umap$CompositionCluster_CC <- NULL
emb.umap$Slide <- NULL
emb.umap <- as.matrix(emb.umap)
colnames(emb.umap) <- c("UMAP1", "UMAP2")

CRC_sp[['umap.ischia12']] <- CreateDimReducObject(embeddings = emb.umap, key = 'umap.ischia12_', assay = 'rctd_tier1')

pdf("8_seurat_ischia_umap_12.pdf")
DimPlot(CRC_sp, reduction = "umap.ischia12", label = FALSE, group.by="cc_12")
dev.off()

# Bar plots
pdf("9_barplot_SampVsorig_12.pdf", height=12, width=20)
dittoBarPlot(CRC_sp, "orig.ident", group.by = "cc_12")
dev.off()

pdf("10_barplot_origVsSamp_12.pdf", height=10, width=20)
dittoBarPlot(CRC_sp, "cc_12", group.by = "orig.ident")
dev.off()

# Cell type co-occurrence
CC4.celltype.cooccur <- spatial.celltype.cooccurence(spatial.object=CRC_sp, deconv.prob.mat=CRC_std, 
                                                     COI="CC4", prob.th= 0.05, 
                                                     Condition=unique(CRC_sp$orig.ident))
pdf("11_celltype_cooccurrence_CC4.pdf")
plot.celltype.cooccurence(CC4.celltype.cooccur)
dev.off()