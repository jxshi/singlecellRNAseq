SubPropPlot <- function(object, groupBy) {
    # Step 1. get data
    plot_data = object@meta.data %>% 
      #dplyr::select({{groupBy}},seurat_clusters) %>% 
      dplyr::select({{groupBy}},subcelltype) %>% 
      dplyr::rename(group = as.name(groupBy), seurat_clusters = as.name("subcelltype"))

    # Step 2. Plot
    figure = ggstatsplot::ggbarstats(data = plot_data, 
                        x = seurat_clusters, y = group,
                        package = 'ggsci',
                        palette = 'category20c_d3',
                        results.subtitle = FALSE,
                        bf.message = FALSE,
                        proportion.test = FALSE,
                        label.args = list(size = 2, 
                                          fill = 'white', 
                                          alpha = 0.85,
                                          family = 'Arial',
                                          fontface = 'bold'),
                                          perc.k = 2,
                                          title = '',
                                          xlab = '',
                                          legend.title = 'Seurat Cluster',
                                          ggtheme = ggpubr::theme_pubclean()) +
      theme(axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(color = 'black', lineend = 'round'),
            legend.position = 'right',
            axis.text.x = element_text(size = 12, color = 'black', family = 'Arial'),
            axis.text.y = element_text(size = 12, color = 'black', family = 'Arial'),
            legend.text = element_text(family = 'Arial', size = 10, color = 'black'),
            legend.title = element_text(family = 'Arial', size = 13, color = 'black')) 
  
    # Step 3. remove some elements
    gginnards::delete_layers(x = figure, match_type = 'GeomText')
}
