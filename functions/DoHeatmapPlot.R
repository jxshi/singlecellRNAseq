DoHeatmapPlot <- function(object, groupBy, features) {
  require(ComplexHeatmap)
  # (1)获取绘图数据
  plot_data = SeuratObject::FetchData(object = object,
                                      vars = c(features, groupBy), 
                                      slot = 'counts') %>% 
    dplyr::mutate(across(.cols = where(is.numeric), .fns = ~ log2(.x + 1))) %>% 
    dplyr::rename(group = as.name(groupBy)) %>% 
    dplyr::arrange(group) %T>% 
    assign(x = 'clusterInfo', value = .$group, envir = .GlobalEnv) %>% 
    dplyr::select(-group) %>% 
    t()
  
  # (2)配色方案：
  require(BuenColors)
  col = jdb_color_maps[1:length(unique(clusterInfo))]
  names(col) = as.character(unique(clusterInfo))
  
  # (3)列(聚类)注释：
  # HeatmapAnnotation函数增加聚类注释，如条形图，点图，折线图，箱线图，密度图：
  top_anno = HeatmapAnnotation(
    cluster = anno_block(gp = gpar(fill = col), #设置填充色；
                         labels = as.character(unique(clusterInfo)), 
                         labels_gp = gpar(cex = 0.5, 
                                          col = 'white',
                                          family = 'Arial'))) #设置字体；
  
  # (4)行注释：
  # rowAnnotation和anno_mark突出重点基因：
  plot_data = as.data.frame(plot_data) 
  gene_pos = which(rownames(plot_data) %in% features) 
  #which和%in%联用返回目标向量位置；
  plot_data = as.matrix(plot_data)
  
  row_anno = rowAnnotation(mark_gene = anno_mark(at = gene_pos, 
                                                 labels = features))
  
  # (5)绘图：
  Heatmap(matrix = plot_data,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_column_names = FALSE,
          show_row_names = FALSE,
          show_heatmap_legend = TRUE,
          column_split = clusterInfo,
          top_annotation = top_anno, #在热图上边增加注释；
          column_title = NULL,
          right_annotation = row_anno,
          use_raster = FALSE,
          heatmap_legend_param = list(
            title = 'log2(count+1)',
            title_position = 'leftcenter-rot'))
}
