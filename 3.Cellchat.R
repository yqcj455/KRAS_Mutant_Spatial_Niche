#### adata转Seurat ####
library(SeuratDisk)
library(Seurat)
name="G:/E-disk/ZJH_KRAS/Workdir/2.Scanpy/ZJH_KRAS_all_"
Convert(paste0(name,"afternewtype_raw.h5ad"), dest = "h5seurat", overwrite = F)
sc_new <- LoadH5Seurat(paste0(name,"afternewtype_raw.h5seurat"),meta.data = F,misc=F)
meta=read.table(paste0(name,"afternewtype_raw_obs_infor.csv"),header = T,sep = ",")

for (col in colnames(meta)) {
  sc_new@meta.data[[col]] <- meta[[col]]
}

scRNA_harmony=sc_new
#随机抽取20000个细胞进行下一步分析
set.seed(23456)  # 为了可重复性，可以设置随机种子
random_cells <- sample(rownames(meta), 20000)
meta=meta[random_cells,]
scRNA_harmony <- subset(scRNA_harmony, cells = meta$X)
table(scRNA_harmony@meta.data[["Maintype"]])
celltype=unique(scRNA_harmony@meta.data[["Maintype"]])
scRNA_harmony <- subset(scRNA_harmony, subset = Maintype %in% c('B_cells', 'Endothelial_cells', 'Epithelial_cells',
                                                                'Fibroblasts', 'Macrophages', 'Monocytes', 
                                                                'Plasma_cells', 'T_CD4', 'T_CD8'))

# 列出当前环境中的所有对???
all_objects <- ls()
# 移除除了 "data", "pd", "fd" 以外的所有对???
rm(list = all_objects[!all_objects %in% c("scRNA_harmony")])
# 调用垃圾收集器，尝试回收未使用的内存
gc()
setwd("G:/E-disk/ZJH_KRAS/Workdir/3.Cellchat")
library(CellChat)
library(tidyverse)
library(ggalluvial)
library(Seurat)
library(data.table)
library(ggsci)
library(SeuratDisk)

#LoadLoom路径不能有中文，不然报错
##直接h5ad转换成h5seurat 
##overwrite参数：覆盖源文件
####可选CellChatDB.human, CellChatDB.mouse
CellChatDB <- CellChatDB.human
##下一步不出图的时候运??? dev.new()
showDatabaseCategory(CellChatDB)

##查看CellchatDB信息
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)

########在CellChat中，我们还可以先择特定的信息描述细胞间的相互作用???
##可以理解为从特定的侧面来刻画细胞间相互作用，比用一个大的配体库又精细了许多???
##查看可以选择的侧???
unique(CellChatDB$interaction$annotation)

#???10万个细胞，随机抽???1/5
sc.con=subset(scRNA_harmony, subset = KRAS =="Wildtype")
sc.exp=subset(scRNA_harmony, subset = KRAS =="Mutant")

#创建Cellchat对象
cellchat.sccon <- createCellChat(object =sc.con@assays$RNA@data, meta =sc.con@meta.data,  group.by ="Maintype")
cellchat.scexp <- createCellChat(object =sc.exp@assays$RNA@data, meta =sc.exp@meta.data,  group.by ="Maintype")
#save(cellchat.sccon ,cellchat.scexp,file = "YS_exp_con_cellchat.rda")
name="ZC_TAM_"

dir.create("Cellchat_Newtype")
setwd("Cellchat_Newtype/")
plan("multisession", workers = 8)
options(future.globals.maxSize = 15*1024^3)
#### cellchat对照??? ####
cellchat=cellchat.sccon 
cellchat@DB  <- subsetDB(CellChatDB, search = c("Secreted Signaling", "ECM-Receptor","Cell-Cell Contact")) # use Secreted Signaling
cellchat <- subsetData(cellchat)
#下面几步时间较长,5000个细胞大???30分钟
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc_con = cellchat
#### cellchat实验???#####
cellchat=cellchat.scexp
cellchat@DB  <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact")) # use Secreted Signaling
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc_exp = cellchat
#### 合并cellchat #####

#这样更新就成功了
# object.list=list(con = object.list[[2]],exp = object.list[[1]])
object.list=list(Wildtype = cc_con,Mutant = cc_exp) #实验组在后，对照组在前
cellchat=mergeCellChat(object.list,cell.prefix = T,add.names = names(object.list))
save(cellchat,object.list,file = paste0(name,"all_cellchat_merged_newtype.rdata")) #这个更重???
##可视???
##所有细胞群总体观：通讯数量与强度对???
pdf(paste0(name,"compareInteractions_count.pdf"))
compareInteractions(cellchat,show.legend = F,group = c(1,2),measure = "count")
dev.off()
pdf(paste0(name,"compareInteractions_weight.pdf"))
compareInteractions(cellchat,show.legend = F,group = c(1,2),measure = "weight")
dev.off()##第一个图展示通讯数量之间的差异，第二个图展示通讯强度之间的差?????? 

#### 数量与强度差异网络图 ####
pdf(paste0(name,"netVisual_diffInteraction_count.pdf"))
netVisual_diffInteraction(cellchat, weight.scale = TRUE)
dev.off()

pdf(paste0(name,"netVisual_diffInteraction_weight.pdf"))
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
dev.off()
##红色是case相对于_con上调的，蓝色是下调的

#### 细胞交互热图 ####
#红色或蓝色表示第二个数据集中与第一个数据集相比增加或[减少]信号
pdf(paste0(name,"netVisual_heatmap_count.pdf"))
netVisual_heatmap(cellchat)
dev.off()

pdf(paste0(name,"netVisual_heatmap_weight.pdf"))
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

library(ggplot2)
#### 保守和特异性信号通路的识别与可视化柱状图 ####
pdf(paste0(name,"rankNet_comparison_stacked_polished.pdf"),height = 5,width = 14)
rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE,do.flip = F)
dev.off()

pdf(paste0(name,"rankNet_comparison.pdf"),height = 14,width = 14)
rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE)
dev.off()
##左图最下面多个信号通路是case组独有的

#### 细胞互作数量对比网络??? ####
weight.max=getMaxWeight(object.list,attribute = c("idents","count"))
pdf(paste0(name,"netVisual_circle_exp_count.pdf"))
netVisual_circle(object.list[[1]]@net$count, weight.scale = TRUE, label.edge = FALSE,
                 edge.weight.max = weight.max[2], edge.width.max = 12, title.name = "")
dev.off()

pdf(paste0(name,"netVisual_circle_con_count.pdf"))
netVisual_circle(object.list[[2]]@net$count, weight.scale = TRUE, label.edge = FALSE,
                 edge.weight.max = weight.max[2], edge.width.max = 12, title.name = "")
dev.off()

weight.max=getMaxWeight(object.list,attribute = c("idents","weight"))
pdf(paste0(name,"netVisual_circle_exp_weight.pdf"))
netVisual_circle(object.list[[1]]@net$weight, weight.scale = TRUE, label.edge = FALSE,
                 edge.weight.max = weight.max[2], edge.width.max = 12, title.name = "")
dev.off()

pdf(paste0(name,"netVisual_circle_con_weight.pdf"))
netVisual_circle(object.list[[2]]@net$weight, weight.scale = TRUE, label.edge = FALSE,
                 edge.weight.max = weight.max[2], edge.width.max = 12, title.name = "")
dev.off()

#### 进入和出去的通路Heatmap模式
#修改行列名这个可以直接反???
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 20)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 20)
pdf(paste0(name,"SignalingPattern_outgoing.pdf"),height = 10,width = 16)
ht1 + ht2
dev.off()

ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 20, color.heatmap = "GnBu")
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 20, color.heatmap = "GnBu")
pdf(paste0(name,"SignalingPattern_incoming.pdf"),height = 10,width = 16)
ht3 + ht4
dev.off()

ht5 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 20, color.heatmap = "OrRd")
ht6 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 20, color.heatmap = "OrRd")
pdf(paste0(name,"SignalingPattern_overall.pdf"),height = 10,width = 16)
ht5 + ht6
dev.off()
#### 配体受体对调节的通信概率气泡??? ####
table(cellchat@meta[["Newtype"]])
commProb <- subsetCommunication(cellchat, 
                                sources.use  = c(12,13,20), 
                                targets.use  = c(12,13,20),
                                thresh = 0.001)
commProb_treat=as.data.frame(commProb$Mutant)
# 按概率降序排序，并提取前30行（含去重处理）
topPathways <- commProb_treat[order(-commProb_treat$prob), ] %>% 
  distinct(pathway_name, .keep_all = TRUE) %>%  # 去重（按通路名）
  head(n = 30)  # 取Top30 
topLRs=as.data.frame(topPathways$interaction_name)
colnames(topLRs)="interaction_name"
#sources.use是来源既发出配体的细胞，targets.use 是受体细???,comparison是组间比???
pdf(paste0(name,"LR_bubble_Epi01_S100A8_CTHRC1_simple.pdf"),height = 20,width = 20)
netVisual_bubble(cellchat, sources.use = c(12,13,20), targets.use = c(12,13,20),  
                 pairLR.use=topLRs,angle.x = 90,thresh = 0.001, remove.isolate = TRUE,
                 comparison = c(1, 2))
dev.off()

#### 提取差异最显著的30个受配体对 ####
wildtype_data <- as.data.frame(commProb$Wildtype)
mutant_data <- as.data.frame(commProb$Mutant)

# 合并两个数据框，并添加分组标识列
combined_data <- rbind(
  cbind(wildtype_data, Group = "Wildtype"),
  cbind(mutant_data, Group = "Mutant")
)

# 计算每个受配体对在两个分组之间的差异
diff_data <- combined_data %>%
  group_by(interaction_name) %>%
  summarize(
    wildtype_prob = mean(prob[Group == "Wildtype"], na.rm = TRUE),
    mutant_prob = mean(prob[Group == "Mutant"], na.rm = TRUE),
    diff_abs = abs(wildtype_prob - mutant_prob)
  )

# 按差异降序排序，并提取前30行
topPathways <- diff_data %>%
  arrange(desc(diff_abs)) %>%
  head(n = 30)

# 提取 Top LRs
topLRs <- as.data.frame(topPathways$interaction_name)
colnames(topLRs) <- "interaction_name"

pdf(paste0(name,"LR_bubble_Epi01_S100A8_CTHRC1_diff_top30.pdf"),height = 10,width = 10)
netVisual_bubble(cellchat, sources.use = c(12,13,20), targets.use = c(12,13,20),  
                 pairLR.use=topLRs,angle.x = 90,thresh = 0.001, remove.isolate = TRUE,
                 comparison = c(1, 2))
dev.off()


subset_data <- subsetCommunication(cellchat, sources.use = c(2,9), targets.use = c(2,9))
# 保存 WT_Liver 数据??? CSV
subset_data2=rbind(subset_data$WT_Liver,subset_data$Parkin_KO_Liver)
library(writexl)
write_xlsx(subset_data2,paste0(name,"LR_bubble.xlsx"))

pdf(paste0(name,"Collagen_bubble2.pdf"),height = 10,width = 16)
netVisual_bubble(cellchat, sources.use = c(1,2,3,5), targets.use = c(8,10,12,13,15),  
                 comparison = c(1, 2), angle.x = 90,signaling = "COLLAGEN")
dev.off()

#此外，我们可以在一个数据集中识别与另一个数据集相比上升（增加）和下降调节（减少）信号配体受体对
#。这可以通过指定max.dataset和min.dataset在函数netVisual_bubble中完成。信号增加意味着这些信号在一个数据集中与其他数据集相比具有更高的通信概率（强度）???
gg1 <- netVisual_bubble(cellchat, sources.use = c(1:10), targets.use = c(4),  
                        comparison = c(1, 2), max.dataset = 2, 
                        title.name = "Increased signaling in _exp", angle.x = 90, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat,sources.use = c(1:10), targets.use = c(4), 
                        comparison = c(1, 2), max.dataset = 1, 
                        title.name = "Decreased signaling in _con", angle.x = 90, remove.isolate = T)
#> Comparing communications on a merged object
pdf(paste0(name,"diff_LR_bubbleplot_All_to_Cd8.pdf"))
gg1 + gg2
dev.off()
#气泡图中显示的配体受体对可以通过signaling.LSIncreased = gg1$data 访问???
#### 挑选目标的受配体对可视??? #### 
# 定义一个阳性数据集，即与其他数据集相比具有正向倍数变化（Fold Change）的数据集
pos.dataset = "Mutant" 
# 定义一个字符名称，用于存储差异表达分析的结果
features.name = pos.dataset 

# 执行差异表达分析，时间很长
cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "datasets", # 按照数据集分组进行分析
                                       pos.dataset = pos.dataset, # 阳性数据集定义为 "LS"
                                       features.name = features.name, # 结果以 "LS" 为名称存储
                                       only.pos = FALSE, # 同时分析正向和负向的表达变化
                                       thresh.pc = 0.1, # 基因需要在至少10%的细胞中表达
                                       thresh.fc = 0.1, # 倍数变化阈值设定为 0.1
                                       thresh.p = 1) # 显著性水平阈值设为 1（较宽松，不严格限制）

#> 使用来自合并 CellChat 对象的联合细胞标签进行分析

# 将差异表达分析的结果映射到推断的细胞间通讯网络中，以便更方便地管理和筛选感兴趣的配体-受体对
net <- netMappingDEG(cellchat, features.name = features.name)

# 提取在数据集 LS 中配体上调的配体-受体对
net.up <- subsetCommunication(cellchat, 
                              net = net, # 使用映射后的网络结果
                              datasets = "Mutant", # 目标数据集为 LS
                              ligand.logFC = 1, # 配体的倍数变化（Log Fold Change）需要大于等于 0.2
                              receptor.logFC = NULL) # 对受体的倍数变化不设限制

# 提取在数据集 NL 中受体和配体均上调（即在 LS 中下调）的配体-受体对
net.down <- subsetCommunication(cellchat, 
                                net = net, # 使用映射后的网络结果
                                datasets = "Wildtype", # 目标数据集为 NL
                                ligand.logFC = -0.3, # 配体的倍数变化需要小于等于 -0.1
                                receptor.logFC = -0.3) # 受体的倍数变化需要小于等于 -0.1#由于信号基因在多亚单位中可能很复杂，我们可以使用net.upnet.down进一步的来获得单个信号基??????

#由于信号基因在多亚单位中可能很复杂，我们可以使用net.upnet.down进一步的来获得单个信号基因。
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
#然后，我们使用气泡图或和弦图可视化上调和向下调的信号配体??????
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, 
                        sources.use = c(1,12:14,20:21,30), targets.use = c(1:5,9,12:13,20:21,30), 
                        comparison = c(1, 2),  angle.x = 90, 
                        remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, 
                        sources.use = c(1:10), targets.use = c(4), 
                        comparison = c(1, 2),  angle.x = 90, 
                        remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pdf(paste0(name,"select_LR_bubbleplot.pdf"))
gg1 + gg2
dev.off()

#### 使用和弦图可视化上调和下调的信号配体??? ####
#par(mfrow = c(1,2), xpd=TRUE) #xpd=TRUE可以允许绘图对象在绘图设备的边界外部分显???
pdf(paste0(name,"netup_chordplot.pdf"))
netVisual_chord_gene(object.list[[2]], sources.use = c(1:10), targets.use = c(4:10),  
                     slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, 
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
dev.off()
pdf(paste0(name,"netdown_chordplot.pdf"))
netVisual_chord_gene(object.list[[1]], sources.use = c(1:10), targets.use = c(4:10), 
                     slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, 
                     title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()

#### 感兴趣的通路可视???
pathways.show <- c("TGFb") #这里只能写一???
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # _con the edge weights across different datasets
pdf(paste0(name,pathways.show,"_netVisual_diffInteraction_count.pdf"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
pdf(paste0(name,pathways.show,"_Heatmap_diffInteraction_count.pdf"),width = 12,height = 6)
#ht_gap = unit(0.5, "cm")则是指定了图形中水平文本之间的间隔为0.5厘米。gap = unit(0.5, "cm"))
dev.off()

# Chord diagram
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
pdf(paste0(name,pathways.show,"_Chord_diffInteraction_count.pdf"),width = 12,height = 6)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

#比较不同数据集之间的信号基因表达分布
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NL", "LS")) # set factor level
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T)
save(cellchat,file = paste0(name,"all_cellchat_afterdiff.rdata")) #这个更重???

####