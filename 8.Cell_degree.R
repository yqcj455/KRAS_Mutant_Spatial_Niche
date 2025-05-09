library(Seurat)
library(dplyr)
library(dbscan)
library(ggplot2)
library(tidyverse)

# 读取数据
ZJH_KRAS_combined_seurat = readRDS('G:/E-disk/ZJH_KRAS/Workdir/6.Cluster/output/ZJH_KRAS_combined_seurat.rds')
samples = unique(ZJH_KRAS_combined_seurat@meta.data[["sample"]])
CRC_cell2loc = read.csv('G:/E-disk/ZJH_KRAS/Workdir/4.Cell2location/ZJH_KRAS_spatial_st_cell2location_res.csv', header = T, row.names = 1, check.names = F)

# 对整个数据框进行 0-1 标准化
CRC_std <- as.data.frame(apply(CRC_cell2loc, 2, function(x) {
  (x - min(x)) / (max(x) - min(x))
}))

# 获取列名
rownames <- rownames(CRC_std)

# 循环修改列名
for (i in seq_along(samples)) {
  # 构造匹配的正则表达式
  pattern <- paste0(samples[i], "_")
  
  # 匹配列名并修改：添加前缀，删除 "_i"
  rownames <- ifelse(
    grepl(pattern, rownames),
    paste0(str_replace(rownames, pattern, ""), "_", i),
    rownames
  )
}

rownames(CRC_std) <- rownames
#sample_name=samples[1]
# 遍历每个样本
for (sample_name in samples) {
  CRC_sp=ZJH_KRAS_combined_seurat[,rownames]
  # 根据当前样本选择数据
  CRC_sp = subset(CRC_sp, subset = sample  == sample_name)
  # 筛选出当前样本的标准化矩阵
  CRC_std_sample = CRC_std[ colnames(CRC_sp),]
  # 选择细胞类型
  decon_mtrx = CRC_std_sample
  cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "max")]
  decon_df <- decon_mtrx %>%
    data.frame(check.names = F) %>%
    tibble::rownames_to_column("barcodes")
  
  CRC_sp@meta.data <- CRC_sp@meta.data %>%
    tibble::rownames_to_column("barcodes") %>%
    dplyr::left_join(decon_df, by = "barcodes") %>%
    tibble::column_to_rownames("barcodes")
  
  # 获取当前样本的组织图像路径
  img_path <- paste0('G:/E-disk/ZJH_KRAS/Spatial_Datadir/', sample_name, '/spatial/tissue_lowres_image.png')
  img <- png::readPNG(img_path)  # 读取低精度图片
  
  img_grob <- grid::rasterGrob(img, interpolate = FALSE, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
  
  # plot dot
  slice <- sample_name
  metadata_ds <- data.frame(CRC_sp@meta.data)
  colnames(metadata_ds) <- colnames(CRC_sp@meta.data)
  cell_types_interest <- cell_types_all
  
  metadata_ds <- metadata_ds %>% tibble::rownames_to_column("barcodeID") %>%
    dplyr::mutate(rsum = base::rowSums(.[, cell_types_interest, drop = FALSE])) %>%
    dplyr::filter(rsum != 0) %>%
    dplyr::select("barcodeID") %>%
    dplyr::left_join(metadata_ds %>%
                       tibble::rownames_to_column("barcodeID"), by = "barcodeID") %>%
    tibble::column_to_rownames("barcodeID")
  
  spatial_coord <- data.frame(CRC_sp@images[[slice]]@coordinates) %>%
    tibble::rownames_to_column("barcodeID") %>%
    dplyr::mutate(imagerow_scaled = imagerow * CRC_sp@images[[slice]]@scale.factors$lowres, 
                  imagecol_scaled = imagecol * CRC_sp@images[[slice]]@scale.factors$lowres) %>%
    dplyr::inner_join(metadata_ds %>%
                        tibble::rownames_to_column("barcodeID"), by = "barcodeID")
  
  st.object = CRC_sp
  st.object$pixel_x = spatial_coord$imagerow_scaled
  st.object$pixel_y = spatial_coord$imagecol_scaled
  
  xys = setNames(st.object@meta.data[, c("pixel_x", "pixel_y", "sample")], c("x", "y", "sample"))
  
  platforms = 'Visium'
  
  spotnames <- rownames(xys)
  names(spotnames) <- c(1:nrow(xys)) %>% paste0()
  
  sdist = 200
  nNeighbours = 6
  maxdist = 200
  
  knn_spatial <- dbscan::kNN(x = xys[, c("x", "y")] %>% as.matrix(), k = nNeighbours)
  knn_spatial.norm <- data.frame(from = rep(1:nrow(knn_spatial$id), nNeighbours),
                                 to = as.vector(knn_spatial$id),
                                 weight = 1/(1 + as.vector(knn_spatial$dist)),
                                 distance = as.vector(knn_spatial$dist))
  
  minK = 4
  spatnet <- knn_spatial.norm
  spatnet$from <- spotnames[spatnet$from]
  spatnet$to <- spotnames[spatnet$to]
  spatnet <- spatnet %>% group_by(from) %>% mutate(rnk = rank(distance)) %>% ungroup()
  spatnet = subset(spatnet, distance <= maxdist | rnk <= minK)
  
  spatnet <- cbind(spatnet, setNames(xys[spatnet$from, 1:2], paste0("start_", c("x", "y"))))
  spatnet <- cbind(spatnet, setNames(xys[spatnet$to, 1:2], paste0("end_", c("x", "y"))))
  
  # 创建输出文件夹
  output_folder <- paste0("ZJH_KRAS_", sample_name, "_degree")
  dir.create(output_folder, showWarnings = FALSE)
  #celltype="Epi_01"
  # 遍历每个细胞类型
  for (celltype in cell_types_interest) {
    threshold = 0.1
    data = as.data.frame(cbind(Barcode = rownames(decon_mtrx), decon_mtrx[, celltype]))
    rownames(data) = data$Barcode
    colnames(data)[2] = celltype
    
    # 提取含有特定细胞类型的spots
    celltype_exist = data[which(data[celltype] > threshold), ]
    
    # 如果没有相关的细胞类型数据，跳过
    if (nrow(celltype_exist) == 0) {
      next
    }
    
    # 创建空的degree矩阵
    degree = matrix(nrow = dim(celltype_exist)[1], ncol = 1)
    rownames(degree) = celltype_exist$Barcode
    colnames(degree) = 'degree'
    
    for (i in rownames(celltype_exist)) {
      su = spatnet[which(spatnet$from == i), ]
      degree[i, 1] = sum(data[su$to, ][celltype] > threshold)
    }
    
    degree = na.omit(degree)
    object = st.object[, rownames(degree)]
    Idents(object) = degree
    
    # 绘图
    Idents(object) = factor(Idents(object), levels = c(0, 1, 2, 3, 4, 5, 6))
    spatial_coord_filtered = spatial_coord[which(spatial_coord$barcodeID %in% rownames(degree)), ]
    spatial_coord_filtered$degree = Idents(object)
    
    pdf_file_path <- file.path(output_folder, paste0("ZJH_KRAS_", sample_name, "_", celltype, "_", threshold, ".pdf"))
    pdf(pdf_file_path,width = 7.5,height = 5)
    p=ggplot2::ggplot() +
      ggplot2::annotation_custom(grob = img_grob, xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
      geom_segment(data = spatnet, aes(x = start_y, y = start_x, xend = end_y, yend = end_x), color = 'white') +
      ggplot2::geom_point(data = spatial_coord_filtered, ggplot2::aes(x = imagecol_scaled, y = imagerow_scaled, color = degree), size = 1.2) +
      ggplot2::ylim(nrow(img), 0) +
      ggplot2::xlim(0, ncol(img)) +
      cowplot::theme_half_open(11, rel_small = 1) +
      ggplot2::theme_void() +
      ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
      scale_color_manual(values = c('#F6FBFC', '#E5F4F8', '#CBEBE5', '#97D7C8', '#67C1A4', '#42AD76', '#258B44')) +
      guides(color = guide_legend(override.aes = list(size = 5))) +
      ggplot2::scale_y_reverse() +
      ggtitle(celltype) +
      theme(plot.title = element_text(size = 20, hjust = 0.5)) +
      guides(size = "none")
    print(p,newpages=F)
    dev.off()
  }
}