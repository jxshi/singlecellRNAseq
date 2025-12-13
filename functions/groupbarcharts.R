groupbarcharts <- function(object, groupBy) {
    # Step 1. get data
    plot_data = object@meta.data %>% 
      dplyr::select({{groupBy}}, celltype) %>% 
      dplyr::rename(group = as.name(groupBy))
    
    # Combine "group" and "celltype" columns into a single column and count occurrences
    result <- plot_data %>%
      group_by(group, celltype) %>%
      summarize(count = n()) %>%
      ungroup() %>%
      group_by(group) %>%
      mutate(total_count = sum(count)) %>%
      mutate(percentage = (count / total_count) * 100)

    # Rename the columns for clarity
    colnames(result) <- c("group", "celltype", "count", "total_count", "percentage")
    
    # write the result data frame
    write.table(result, "celltype_percentage_in_each_group.tsv", quote = F, sep="\t", row.names=F)
    # Step 2. Plot
    p <-  ggplot(data=result, aes(x = celltype, y= percentage, fill = group)) + 
              geom_bar(stat = "identity", position = "dodge") +
  
              # Customize the plot
              labs(
                title = "Percentage of Cell Types in Each Group",
                x = "Group",
                y = "Percentage"
              ) +
      
              theme_minimal() +
     # Remove background and grid lines
            theme(
    		panel.background = element_blank(),
    		panel.grid.major = element_blank(),
    		panel.grid.minor = element_blank(),
    		axis.line = element_line(color = "black")
  		) +
	    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
    p <- p + ggsci::scale_fill_nejm() + coord_flip()
    p
}
