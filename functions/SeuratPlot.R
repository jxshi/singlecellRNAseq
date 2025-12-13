SeuratPlot <- function(object, groupBy, reduction){
  # groupBy = 'celltype'
  # object = sce
  
  # (1)获取非线性降维坐标，可以选择tsne或者umap，前提是存在哦：
  plot_data = object@reductions[[reduction]]@cell.embeddings %>% 
    as.data.frame() %>% 
    dplyr::mutate(Barcode = rownames(.)) %>% 
    inner_join(object@meta.data %>% rownames_to_column('Barcode'), by = 'Barcode') %>% 
    column_to_rownames('Barcode') %>% 
    dplyr::select(1,2, {{groupBy}}) %>% 
    dplyr::rename(group = as.name(groupBy))
  
  # (2)生成聚类中心坐标
  centroids = aggregate(as.matrix(plot_data[,c(1,2)]) ~ group,
                        data = plot_data,
                        FUN = mean)

  # (3)绘图
  if(reduction == 'tsne') {
    x = 'tSNE_1'; y = 'tSNE_2'
  } else {
    x = 'UMAP_1'; y = 'UMAP_2'
  }
  ggplot(data = plot_data, mapping = aes(x = !!as.name(x), y = !!as.name(y))) +
    geom_point(mapping = aes(fill = group),
               color = 'white', shape = 21, stroke = 0.5, size = 3) +
    geom_point(data = centroids,
               mapping = aes(x = !!as.name(x), y = !!as.name(y),
                             fill = group),
               color = 'white', shape = 21, stroke = 0.5, size = 3) +
    theme_bw() +
    ggrepel::geom_label_repel(data = centroids,
                              mapping = aes(label = group,
                                            color = group),
                              size = 7, fill = alpha('white', 0.6),
                              fontface = 'bold',
                              family = 'Arial',
                              show.legend = FALSE) +
    scale_fill_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.65)) +
    scale_color_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.65)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(family = 'Arial', size = 12, color = 'black'),
          axis.text.y = element_text(family = 'Arial', size = 12, color = 'black'),
          axis.title.x = element_text(family = 'Arial', size = 15, color = 'black'),
          axis.title.y = element_text(family = 'Arial', size = 15, color = 'black'),
          axis.ticks = element_line(color = 'black', lineend = 'round'),
          legend.position = 'bottom',
          legend.text = element_text(family = 'Arial', size = 10, color = 'black'),
          legend.title = element_text(family = 'Arial', size = 13, color = 'black'))
}
