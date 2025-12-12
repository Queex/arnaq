#' Creates a scatter plot of reads for two samples.
#' 
#' The scatter plot compares the reads for genes across two samples, preferably using normalised
#' data. A filter is applied for minimum read depth for a gene, optionally using a different
#' dat matrix (e.g. unnormalised).
#' 
#' @param norm.dat The dataframe containing the reads to base the scatterplot on.
#' @param dat.for.filter A similar dataframe containing the data to use for filtering. By default 
#' the same as `norm.dat`.
#' @param col1name The name of the coolumn containing data for sample 1.
#' @param col2name The name of the coolumn containing data for sample 2.
#' @param title The title to use for the plot.
#' @param log If `TRUE`, use log2 scales for the axes.
#' @param filter.threshold The minimum value to filter on in `dat.for.filter`. The default of 3
#' prevents numerical issues when preparing the plot.
#' @returns a scatterplot.
#' @seealso `\link{plot_all_scatter_pairs}`
#' 
#' @export
plot_reads_scatter <- function(norm.dat, dat.for.filter=norm.dat, col1name, col2name, 
                               title = NULL, log = FALSE, filter.threshold = 3) {
  dat.for.filter <- dat.for.filter[, c(col1name, col2name)]
  dat <- data.frame(norm.dat[, c(col1name, col2name)])
  dat <- dat[rowSums(dat.for.filter >= filter.threshold) == 2, ]

  dat.corr <- round(stats::cor(dat[, 1], dat[, 2]), 3)

  if (is.null(title)) {
    title <- paste(col1name, "v", col2name)
  }

  theplot <- ggplot2::ggplot(dat, ggplot2::aes(x = .data[[col1name]], y = .data[[col2name]])) +
    ggplot2::geom_hex(bins = 50) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(col1name) +
    ggplot2::ylab(col2name) +
    ggplot2::geom_smooth(method = "glm", colour = "red", formula = y ~ x) +
    ggplot2::geom_abline(intercept = 0, slope = 1, colour = "green") +
    ggplot2::labs(caption = paste("Correlation:", dat.corr)) +
    ggplot2::scale_fill_gradientn(colours = c(
      "darkcyan", "cyan", "red3", "orange2", "yellow2",
      "yellow"
    ))
  if (log) {
    theplot <- theplot + ggplot2::scale_x_continuous(trans = "log2") +
      ggplot2::scale_y_continuous(trans = "log2")
  }

  return(theplot)
}

#' Creates a list of scatter plot of reads for two samples.
#' 
#' The scatter plots compare the reads for genes across two samples, preferably using normalised
#' data. A filter is applied for minimum read depth for a gene, optionally using a different
#' dat matrix (e.g. unnormalised).
#' 
#' @param norm.dat The dataframe containing the reads to base the scatterplot on.
#' @param dat.for.filter A similar dataframe containing the data to use for filtering. By default 
#' the same as `norm.dat`.
#' @param scatter.pairs A list of length-2 vectors, each entry in the list giving the names of two 
#' samples to create a scatter plot for.
#' @param log If `TRUE`, use log2 scales for the axes.
#' @param filter.threshold The minimum value to filter on in `dat.for.filter`. The default of 3
#' prevents numerical issues when preparing the plot.
#' @returns a scatterplot.
#' @seealso `\link{plot_all_scatter_pairs}`
#' 
#' @export
plot_all_scatter_pairs <- function(norm.dat, dat.for.filter=norm.dat, scatter.pairs, log = TRUE,
                                   filter.threshold = 3) {
  out <- list()
  for (pair in scatter.pairs) {
    out[[paste(pair[1], "v", pair[2])]] <- plot_reads_scatter(norm.dat, dat.for.filter,
      pair[1], pair[2],
      log = log,
      filter.threshold = filter.threshold
    )
  }
  return(out)
}