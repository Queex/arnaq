#' Create a list of four plots to show similarity between samples.
#' 
#' The plots are: a heatmap based on Euclidean distance, a dendrogram showing hierarchical 
#' clustering based on Euclidean distance, a heatmap based on Spearman's rank correlation of the
#' Euclidean distances, and a dendrogram showing hierarchical clustering based on that Spearman's
#' rank correlation.
#' @param data The table of read counts.
#' @returns A list of the four plots.
#' 
#' @export
plot_all_clustering <- function(data) {
  dists <- dist(t(data))
  spearman.dists <- cor(data, method = "spearman")
  return(list(
    plot_clustering_heatmap(as.matrix(dists), "Distance-based heatmap", reverse.palette = FALSE),
    plot_clustering_dendrogram(hclust(dists), "Distance-based hierarchical clustering"),
    plot_clustering_heatmap(spearman.dists, "Spearman correlation heatmap", reverse.palette = TRUE),
    plot_clustering_dendrogram(
      hclust(as.dist(1 - spearman.dists), method = "ward.D2"),
      "Spearman's rank-based hierarchical clustering"
    )
  ))
}

# Convenience function for producing similarity heatmaps.
plot_clustering_heatmap <- function(mat, main, reverse.palette = FALSE) {

  long.mat <- reshape2::melt(as.matrix(mat)) # long format, symmetric anyway

  if (reverse.palette) {
    rgb.palette <- c("yellow", "firebrick")
  } else {
    rgb.palette <- c("firebrick", "yellow")
  }

  p <- ggplot2::ggplot(long.mat, ggplot2::aes(Var1, Var2)) +
    ggplot2::geom_tile(ggplot2::aes(fill = value), colour = "white") +
    ggplot2::scale_fill_gradient(low = rgb.palette[1], high = rgb.palette[2]) +
    ggplot2::labs(x = "", y = "", title = main) +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.background = ggplot2::element_blank()
    ) +
    ggplot2::coord_fixed(ratio = 1)

  # This will need something to change the order of samples as per dendrogram eventually
  return(p)
}

# Convenience function for producing dendrograms
plot_clustering_dendrogram <- function(hclust.obj, title) {
  dend <- as.dendrogram(hclust.obj)
  the.plot <- ggplot2::ggplot(dendextend::as.ggdend(dend)) + ggplot2::labs(title=title)
  heights <- dendextend::get_branches_heights(dend)
  the.plot <- the.plot + ggplot2::ylim(-2 * min(heights), max(heights))
  return(the.plot)
}