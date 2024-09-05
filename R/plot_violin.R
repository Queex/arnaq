#' Creates two violin plots of metrics, based on scale.
#'
#' This assumes that the `count.metrics` argument is a list with 2 entries: one named `ByBase` for
#' count-related metrics, and one named `ByOther` for proportional metrics. This is in order to
#' keep the latter readable without putting multiple scales on a single plot, which would be
#' difficult to do well in an automated fashion.
#'
#' The violin plots are also separated by a group. `/link{arnaq}` uses the first entry in its
#' `treat.group` argument to produce these, run this function separately and regenerate the output
#' if you want it split by another grouping.
#' @param count.metrics A list with two named entries, `ByBase` and `ByOther`, each containing
#' a dataframe of metrics with samples as columns.
#' @param group The group to split on.
#' @returns A list containing two plots.
#' @seealso `\link{arnaq}`
#' @seealso `\link{arnaq.create.report}`
#'
#' @export
plot_all_metrics_violin <- function(count.metrics, group) {
  cat("Making violin plots\n")
  violin1 <- plot_metrics_violin(
    count.metrics$ByBases,
    sample.metadata[, c(3, group)], c("Category", "Million Bases"),
    "Count Metrics", 1 / 1000000
  )
  violin2 <- plot_metrics_violin(
    count.metrics$ByOther,
    sample.metadata[, c(3, group)], c("Metric", "Value"), "Other Metrics", 1
  )
  return(list(violin1, violin2))
}

#' Creates a violin plot of metrics.
#'
#' This creates a violin plot showing metrics (e.g. produced by Picard) for samples, stratified by
#' group. All metrics in the supplied table will be included in a signle plot.
#'
#' @param metric.table The table of metrics, each column a sample.
#' @param group.data A vector of group membership for the samples.
#' @param labels A vector of length 2, giving the axes labels.
#' @param title The title for the plot.
#' @param yscale What scaling to apply to the data. For example, a `yscale` of 100 can be used if
#' you wnat the y axis units to be percent.
#' @returns A violin plot of the metrics.
#'
#' @export
plot_metrics_violin <- function(metric.table, group.data, labels, title, yscale = 1) {
  metric.table <- metric.table * yscale

  metric.table <- cbind(Sample = rownames(metric.table), metric.table)
  dat <- reshape2::melt(metric.table, id.var = "Sample")
  tmp <- group.data[match(dat$Sample, group.data[, 1]), 2]
  dat <- cbind(dat, Group = tmp)

  theplot <- ggplot2::ggplot(dat, ggplot2::aes(x = variable, y = value, fill = Group)) +
    ggplot2::geom_violin(scale = "width", draw_quantiles = c(0.25, 0.5, 0.75), trim = FALSE) +
    ggplot2::labs(x = labels[1], y = labels[2], title = title) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))

  return(theplot)
}