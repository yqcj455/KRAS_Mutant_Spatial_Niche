library(Seurat)
library(tidyverse)
library(SeuratDisk)

name="G:/E-disk/ZJH_KRAS/Workdir/2.Scanpy/ZJH_KRAS_all_"
Convert(paste0(name,"afternewtype.h5ad"), dest = "h5seurat", overwrite = F)
sc_new <- LoadH5Seurat(paste0(name,"afternewtype.h5seurat"),meta.data = F,misc=F)
meta=read.table(paste0(name,"afternewtype_obs_infor.csv"),header = T,sep = ",")

for (col in colnames(meta)) {
  sc_new@meta.data[[col]] <- meta[[col]]
}
CRC_sc_data=sc_new

#### 1.鍗曠粏鑳炲樊寮傚熀鍥犲鎵? ####
# 鎺у埗涓変釜闃堝€硷細logfc.threshold p_val_adj d
Idents(CRC_sc_data)="Newtype"
plan("multisession", workers = 4)
options(future.globals.maxSize = 15*1024^3)
sc_marker <- FindAllMarkers(CRC_sc_data, logfc.threshold = 0.25, only.pos = T)
sc_marker=FindAllMarkers(CRC_sc_data,logfc.threshold = 0.25,only.pos = T)
sc_marker=sc_marker%>%filter(p_val_adj < 1e-05)
sc_marker$d=sc_marker$pct.1 - sc_marker$pct.2
sc_marker=sc_marker%>%filter(d > 0.2)
sc_marker=sc_marker%>%arrange(cluster,desc(avg_log2FC))
sc_marker=as.data.frame(sc_marker)
save(sc_marker,file="ZJH_KRAS_sc_marker.Rdata")
load("ZJH_KRAS_sc_marker.Rdata")
### 2.鎵归噺杩愯绌洪棿杞綍缁? ####
# 娉ㄦ剰锛氳繖涓枃鐚彧鑳戒笅杞藉埌绌鸿浆鐭╅樀鍜宻pot鍧愭爣锛屾病鏈夊浘鍍忋€傚洜姝や笅闈㈢殑娴佺▼鍜屽崟缁嗚優杞綍缁勪竴妯′竴鏍凤紝
# 鑻ユ槸瀹屾暣鐨勭┖杞暟鎹泦锛屽垯鏈変竴浜涙楠ゆ槸涓嶄竴鏍风殑锛屽叿浣撴祦绋嬪彲浠ュ弬鑰僺eurat瀹樼綉锛宧ttps://satijalab.org/seurat/articles/spatial_vignette.html
# 鑾峰彇鎵€鏈夋牱鏈悕绉?
samples <- list.files("/E-disk/ZJH_KRAS/Spatial_Datadir/")
# 寰幆澶勭悊姣忎釜鏍锋湰
for (sample in samples) {
  # 鍔犺浇鏁版嵁
  CRC_sp = Load10X_Spatial(paste0("G:/E-disk/ZJH_KRAS/Spatial_Datadir/", sample))
  
  # 鏁版嵁棰勫鐞?
  CRC_sp <- SCTransform(CRC_sp, assay = 'Spatial', verbose = FALSE, variable.features.n = 2000)
  CRC_sp <- RunPCA(CRC_sp, assay = "SCT", verbose = FALSE, npcs = 20)
  CRC_sp <- FindNeighbors(CRC_sp, reduction = "pca", dims = 1:20, k.param = 10)
  CRC_sp <- FindClusters(CRC_sp, verbose = FALSE, resolution = 0.3)
  CRC_sp <- RunUMAP(CRC_sp, reduction = "pca", dims = 1:20)
  CRC_sp <- RunTSNE(CRC_sp, reduction = "pca", dims = 1:20)
  
  # 鏍规嵁鑱氱被缁撴灉鐢熸垚鏂扮殑cluster鍚嶇О
  new_cluster_names <- paste0("Niche_", seq_len(length(unique(CRC_sp@meta.data[["seurat_clusters"]]))))
  names(new_cluster_names) <- levels(CRC_sp@meta.data[["seurat_clusters"]])
  
  # 鍒涘缓 'Region' 鍒楀苟璧嬪€兼柊鐨刢luster鍚嶇О
  CRC_sp@meta.data[["Region"]] <- new_cluster_names[as.character(CRC_sp@meta.data[["seurat_clusters"]])]
  colors=c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')
  CRC_sp@meta.data[["Region"]]=factor(CRC_sp@meta.data[["Region"]],levels = paste0("Niche_", seq_len(length(unique(CRC_sp@meta.data[["seurat_clusters"]])))))
  cluster_num=length(unique(CRC_sp@meta.data[["Region"]]))
  color_region=colors[1:cluster_num]
  names(color_region)=levels(CRC_sp@meta.data[["Region"]])
  # 淇濆瓨绌洪棿鍒嗗竷鍥?
  pdf(paste0(sample, "_spatial_regions_plot_newcolor.pdf"), width = 6, height = 5)
  # 缁樺埗绌洪棿鍖哄煙鍥惧苟璋冩暣鐐圭殑澶у皬鍜岄鑹?
  p <- SpatialDimPlot(CRC_sp, group.by = "Region", pt.size.factor = 1.3, cols = color_region) +
    theme(
      legend.title = element_text(hjust = 0.5),
      panel.border = element_rect(colour = "black", fill=NA, size=1),  # 娣诲姞澶栬竟妗?
      axis.title.x = element_text(size = 14),  # 璁剧疆妯潗鏍囩殑鏂囨湰鏍峰紡
      axis.title.y = element_text(size = 14)   # 璁剧疆绾靛潗鏍囩殑鏂囨湰鏍峰紡
    ) +
    labs(
      x = "spatial1",  # 璁剧疆妯潗鏍囨爣绛?
      y = "spatial2"   # 璁剧疆绾靛潗鏍囨爣绛?
    ) +
    guides(col = guide_legend(override.aes = list(shape = 16))) +  # 鍘绘帀鐐圭殑榛戣壊杈规
    scale_color_manual(values = color_region)  # 纭繚浣跨敤鑷畾涔夐鑹?
  # 杈撳嚭鍥惧舰鍒癙DF鏂囦欢
  print(p)
  dev.off()
  
  save(CRC_sp, file = paste0(sample, "_after_regions.Rda"))
}

#### 3.MIA鍒嗘瀽绀轰緥 ####
sample=samples[1]
load(paste0(sample, "_after_regions.Rda"))
# 鎴栬€呮牴鎹煇涓熀鍥犵殑琛ㄨ揪姘村钩缁樺埗绌洪棿鍥撅紝姣斿鍩哄洜"GeneOfInterest"
### 鏁村悎鍧愭爣銆乺egion
#锛堟澶勬妸鍧愭爣绫绘瘮鍗曠粏鑳炶浆褰曠粍鍒嗘瀽涓殑娉ㄩ噴淇℃伅???

### 鐢诲浘鐪嬬湅
library(RColorBrewer)
library(scales)
colors=c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')
CRC_A1_sp = CRC_sp
CRC_A1_sp@meta.data[["Region"]]=factor(CRC_A1_sp@meta.data[["Region"]],levels = paste0("Niche_", seq_len(length(unique(CRC_A1_sp@meta.data[["seurat_clusters"]])))))
cluster_num=length(unique(CRC_A1_sp@meta.data[["Region"]]))
color_region=colors[1:cluster_num]
names(color_region)=levels(CRC_A1_sp@meta.data[["Region"]])

CRC_A1_sp@meta.data["x_coord"]=CRC_A1_sp@images[["slice1"]]@coordinates[["row"]]
CRC_A1_sp@meta.data["y_coord"]=CRC_A1_sp@images[["slice1"]]@coordinates[["col"]]

CRC_A1_sp@meta.data%>%ggplot(aes(x=x_coord,y=y_coord,fill=Region))+
  geom_tile(color="white")+
  scale_fill_manual("Cluster assignments",values = color_region)+
  theme_void()+
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "right",
    legend.direction = "vertical"
  )+
  guides(fill = guide_legend(override.aes = list(size=10)))
ggsave(paste0(sample,"_MIA_Cluster_assign.pdf"),width = 14,height = 10,units = "cm")

### 鎵綬egion鐗瑰紓鍩哄洜
Idents(CRC_A1_sp)="Region"
#闃堝€?0.1-0.2
Region_marker=FindAllMarkers(CRC_A1_sp,logfc.threshold = 0.1,only.pos = T)
Region_marker=Region_marker%>%filter(p_val_adj < 0.01)
Region_marker$d=Region_marker$pct.1 - Region_marker$pct.2
Region_marker=Region_marker%>%filter(d > 0.05)
Region_marker=Region_marker%>%arrange(cluster,desc(avg_log2FC))
Region_marker=as.data.frame(Region_marker)

# 璇存槑
# 1. 涓婅堪鎵句袱涓狣EG鏁版嵁妗嗙殑鏂规硶涓嶅敮涓€锛岄槇鍊间篃涓嶅敮涓€
# 2. 绗簩涓狣EG鏁版嵁妗嗕篃鍙互鏄痗luster鐨刴arker
# 3. MIA鍒嗘瀽妯″紡鍦ㄥ崟缁嗚優鍜岀┖闂磋浆褰曠粍鍦烘櫙閮藉彲浠ュ簲鐢紝绌鸿浆鍦烘櫙鏄湅缁嗚優浜氱兢鐨勫瘜闆嗙▼搴︼紝鍗曠粏鑳炲満鏅槸鍋氱粏鑳炰簹缇ゆ敞閲?

region_specific=Region_marker[,c("cluster","gene")]
colnames(region_specific)[1]="region"

celltype_specific=sc_marker[,c("cluster","gene")]
colnames(celltype_specific)[1]="celltype"

N=length(union(rownames(CRC_sc_data),rownames(CRC_A1_sp)))

library(RColorBrewer)
library(scales)
colors=c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')
color_region=colors[1:cluster_num]
names(color_region)=levels(CRC_A1_sp@meta.data[["Region"]])

source("syMIA.R")
miares=syMIA(region_specific,celltype_specific,N,color_region,plotname=paste0(sample,"_MIAplot"),plotwidth=28,plotheight=22)

#### 4.鎵归噺杩愯MIA ####
library(Seurat)
library(tidyverse)
library(SeuratDisk)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(dplyr)
# 鎸夋牱鏈惊鐜?
samples <- list.files("/E-disk/ZJH_KRAS/Spatial_Datadir/")
# 瀹氫箟棰滆壊鍚戦噺
colors <- c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')
for (sample in samples) {
  # 鍔犺浇鏁版嵁
  load(paste0(sample, "_after_regions.Rda"))
  # 璁剧疆 Seurat 瀵硅薄
  CRC_A1_sp <- CRC_sp
  # 澶勭悊 Region 鍒?
  CRC_A1_sp@meta.data[["Region"]] <- factor(CRC_A1_sp@meta.data[["Region"]], levels = paste0("Niche_", seq_len(length(unique(CRC_A1_sp@meta.data[["seurat_clusters"]])))))
  cluster_num <- length(unique(CRC_A1_sp@meta.data[["Region"]]))
  color_region <- colors[1:cluster_num]
  names(color_region) <- levels(CRC_A1_sp@meta.data[["Region"]])
  
  # 娣诲姞鍧愭爣淇℃伅
  CRC_A1_sp@meta.data["x_coord"] <- CRC_A1_sp@images[["slice1"]]@coordinates[["row"]]
  CRC_A1_sp@meta.data["y_coord"] <- CRC_A1_sp@images[["slice1"]]@coordinates[["col"]]
  
  # 缁樺埗绌洪棿鍥?
  p <- CRC_A1_sp@meta.data %>%
    ggplot(aes(x = x_coord, y = y_coord, fill = Region)) +
    geom_tile(color = "white") +
    scale_fill_manual("Cluster assignments", values = color_region) +
    theme_void() +
    theme(
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.position = "right",
      legend.direction = "vertical"
    ) +
    guides(fill = guide_legend(override.aes = list(size = 10)))
  
  ggsave(paste0(sample, "_MIA_Cluster_assign.pdf"), width = 14, height = 10, units = "cm")
  
  # 鎵? Region 鐗瑰紓鍩哄洜
  Idents(CRC_A1_sp) <- "Region"
  Region_marker <- FindAllMarkers(CRC_A1_sp, logfc.threshold = 0.1, only.pos = TRUE)
  Region_marker <- Region_marker %>% filter(p_val_adj < 0.01)
  Region_marker$d <- Region_marker$pct.1 - Region_marker$pct.2
  Region_marker <- Region_marker %>% filter(d > 0.05)
  Region_marker <- Region_marker %>% arrange(cluster, desc(avg_log2FC))
  Region_marker <- as.data.frame(Region_marker)
  
  # 鍑嗗 MIA 鏁版嵁
  region_specific <- Region_marker[, c("cluster", "gene")]
  colnames(region_specific)[1] <- "region"
  
  celltype_specific <- sc_marker[, c("cluster", "gene")]
  colnames(celltype_specific)[1] <- "celltype"
  
  N <- length(union(rownames(CRC_sc_data), rownames(CRC_A1_sp)))
  
  # 杩涜 MIA 鍒嗘瀽
  source("syMIA.R")
  miares <- syMIA(region_specific, celltype_specific, N, color_region, plotname = paste0(sample, "_MIAplot"), plotwidth = 28, plotheight = 22)
}

#### 5.鍗曠嫭淇敼涓€涓牱鏈? ####
samples <- list.files("/E-disk/ZJH_KRAS/Spatial_Datadir/")
# 寰幆澶勭悊姣忎釜鏍锋湰
sample=samples[4]
# 鍔犺浇鏁版嵁
CRC_sp = Load10X_Spatial(paste0("G:/E-disk/ZJH_KRAS/Spatial_Datadir/", sample))

# 鏁版嵁棰勫鐞?
CRC_sp <- SCTransform(CRC_sp, assay = 'Spatial', verbose = FALSE, variable.features.n = 2000)
CRC_sp <- RunPCA(CRC_sp, assay = "SCT", verbose = FALSE, npcs = 20)
CRC_sp <- FindNeighbors(CRC_sp, reduction = "pca", dims = 1:20, k.param = 10)
CRC_sp <- FindClusters(CRC_sp, verbose = FALSE, resolution = 0.2)
CRC_sp <- RunUMAP(CRC_sp, reduction = "pca", dims = 1:20)
CRC_sp <- RunTSNE(CRC_sp, reduction = "pca", dims = 1:20)

# 鏍规嵁鑱氱被缁撴灉鐢熸垚鏂扮殑cluster鍚嶇О
new_cluster_names <- paste0("Niche_", seq_len(length(unique(CRC_sp@meta.data[["seurat_clusters"]]))))
names(new_cluster_names) <- levels(CRC_sp@meta.data[["seurat_clusters"]])

# 鍒涘缓 'Region' 鍒楀苟璧嬪€兼柊鐨刢luster鍚嶇О
CRC_sp@meta.data[["Region"]] <- new_cluster_names[as.character(CRC_sp@meta.data[["seurat_clusters"]])]
colors=c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')
CRC_sp@meta.data[["Region"]]=factor(CRC_sp@meta.data[["Region"]],levels = paste0("Niche_", seq_len(length(unique(CRC_sp@meta.data[["seurat_clusters"]])))))
cluster_num=length(unique(CRC_sp@meta.data[["Region"]]))
color_region=colors[1:cluster_num]
names(color_region)=levels(CRC_sp@meta.data[["Region"]])
# 淇濆瓨绌洪棿鍒嗗竷鍥?
pdf(paste0(sample, "_spatial_regions_plot_newcolor2.pdf"), width = 6, height = 5)
# 缁樺埗绌洪棿鍖哄煙鍥惧苟璋冩暣鐐圭殑澶у皬鍜岄鑹?
p <- SpatialDimPlot(CRC_sp, group.by = "Region", pt.size.factor = 1.3, cols = color_region) +
  theme(
    legend.title = element_text(hjust = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, size=1),  # 娣诲姞澶栬竟妗?
    axis.title.x = element_text(size = 14),  # 璁剧疆妯潗鏍囩殑鏂囨湰鏍峰紡
    axis.title.y = element_text(size = 14)   # 璁剧疆绾靛潗鏍囩殑鏂囨湰鏍峰紡
  ) +
  labs(
    x = "spatial1",  # 璁剧疆妯潗鏍囨爣绛?
    y = "spatial2"   # 璁剧疆绾靛潗鏍囨爣绛?
  ) +
  guides(col = guide_legend(override.aes = list(shape = 16))) +  # 鍘绘帀鐐圭殑榛戣壊杈规
  scale_color_manual(values = color_region)  # 纭繚浣跨敤鑷畾涔夐鑹?
# 杈撳嚭鍥惧舰鍒癙DF鏂囦欢
print(p)
dev.off()

save(CRC_sp, file = paste0(sample, "_after_regions.Rda"))

library(Seurat)
library(tidyverse)
library(SeuratDisk)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(dplyr)
# 鎸夋牱鏈惊鐜?
samples <- list.files("/E-disk/ZJH_KRAS/Spatial_Datadir/")
# 瀹氫箟棰滆壊鍚戦噺
colors <- c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')

# 鍔犺浇鏁版嵁
load(paste0(sample, "_after_regions.Rda"))
# 璁剧疆 Seurat 瀵硅薄
# 澶勭悊 Region 鍒?
CRC_sp@meta.data[["Region"]] <- factor(CRC_sp@meta.data[["Region"]], levels = paste0("Niche_", seq_len(length(unique(CRC_sp@meta.data[["seurat_clusters"]])))))
cluster_num <- length(unique(CRC_sp@meta.data[["Region"]]))
color_region <- colors[1:cluster_num]
names(color_region) <- levels(CRC_sp@meta.data[["Region"]])

# 娣诲姞鍧愭爣淇℃伅
CRC_sp@meta.data["x_coord"] <- CRC_sp@images[["slice1"]]@coordinates[["row"]]
CRC_sp@meta.data["y_coord"] <- CRC_sp@images[["slice1"]]@coordinates[["col"]]

# 缁樺埗绌洪棿鍥?
p <- CRC_sp@meta.data %>%
  ggplot(aes(x = x_coord, y = y_coord, fill = Region)) +
  geom_tile(color = "white") +
  scale_fill_manual("Cluster assignments", values = color_region) +
  theme_void() +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "right",
    legend.direction = "vertical"
  ) +
  guides(fill = guide_legend(override.aes = list(size = 10)))

ggsave(paste0(sample, "_MIA_Cluster_assign.pdf"), width = 14, height = 10, units = "cm")

# 鎵? Region 鐗瑰紓鍩哄洜
Idents(CRC_sp) <- "Region"
Region_marker <- FindAllMarkers(CRC_sp, logfc.threshold = 0.1, only.pos = TRUE)
Region_marker <- Region_marker %>% filter(p_val_adj < 0.1)
Region_marker$d <- Region_marker$pct.1 - Region_marker$pct.2
Region_marker <- Region_marker %>% filter(d > 0.05)
Region_marker <- Region_marker %>% arrange(cluster, desc(avg_log2FC))
Region_marker <- as.data.frame(Region_marker)

# 鍑嗗 MIA 鏁版嵁
region_specific <- Region_marker[, c("cluster", "gene")]
colnames(region_specific)[1] <- "region"

celltype_specific <- sc_marker[, c("cluster", "gene")]
colnames(celltype_specific)[1] <- "celltype"

N <- length(union(rownames(CRC_sc_data), rownames(CRC_sp)))

# 杩涜 MIA 鍒嗘瀽
source("syMIA.R")
miares <- syMIA(region_specific, celltype_specific, N, color_region, plotname = paste0(sample, "_MIAplot"), plotwidth = 28, plotheight = 22)
