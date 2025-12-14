groupbarchartssubcluster <- function(object, groupBy) {
    if (missing(groupBy)) {
        stop("`groupBy` must be provided and refer to a metadata column.")
    }
    if (!"subcelltype" %in% colnames(object@meta.data)) {
        stop("Metadata column `subcelltype` is required for groupbarchartssubcluster().")
    }
    if (!as.character(substitute(groupBy)) %in% colnames(object@meta.data)) {
        stop("`groupBy` must exist in object@meta.data.")
    }

    plot_data <- object@meta.data %>%
      dplyr::select({{groupBy}}, subcelltype) %>%
      dplyr::rename(group = as.name(groupBy)) %>%
      tidyr::drop_na()

    result <- plot_data %>%
      dplyr::group_by(group, subcelltype) %>%
      dplyr::summarize(count = dplyr::n(), .groups = "drop_last") %>%
      dplyr::mutate(total_count = sum(count), percentage = (count / total_count) * 100) %>%
      dplyr::ungroup()

    colnames(result) <- c("group", "subcelltype", "count", "total_count", "percentage")

    write.table(result, "subcelltype_percentage_in_each_group.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

    ggplot(data = result, aes(x = subcelltype, y = percentage, fill = group)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(
        title = "Percentage of Cell Types in Each Group",
        x = "Group",
        y = "Percentage"
      ) +
      theme_minimal() +
      theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      ggsci::scale_fill_nejm()
}
