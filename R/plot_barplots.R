#' Create read count barplots
#' 
#' From a read assignment table, produces a list of plots for the samples in the table.
#' @param counts The read count table, samples as columns.
#' @param max.samples.per.page The maximum number of bars a single plot should have. When this is
#' exceeded, there will be multiple plots produced.
#' @param proportional If `TRUE`, the plots will show the proportion of each assignment category
#' for that sample instead of the absolute number.
#' @returns A list of plots.
#' 
#' @export
plot_counts_barplot <- function(counts, title="Read assignment by category", 
                                max.samples.per.page = 40, proportional = FALSE) {
  if (proportional) {
    fill.position <- TRUE
    ylbl <- "Proportion of reads"
  } else {
    fill.position <- FALSE
    counts <- counts / 10^6
    ylbl <- "Million reads"
  }

  return(plot_barplot_helper_multipage(
    counts, max.samples.per.page,
    fill.position, title, "Sample", ylbl
  ))
}

# Handles the pagination of barplots.
plot_barplot_helper_multipage <- function(counts, max.samples.per.page,
                                          fill.position, title, xlbl, ylbl) {
  counts <- counts[rev(seq_len(nrow(counts))), ]
  n.samples <- ncol(counts)
  pages <- ceiling(n.samples / max.samples.per.page)
  max.samples.per.page <- ceiling(n.samples / pages)
  min.n <- 1

  # Set up y axis range
  y.lim <- c(0, max(colSums(counts)))

  if (fill.position) {
    position <- "fill"
  } else {
    position <- "stack"
  }

  out <- list()
  i <- 1
  while (min.n <= n.samples) {
    max.n <- min.n + max.samples.per.page - 1
    tmp.counts <- counts[, min.n:min(max.n, n.samples)]
    out[[i]] <- plot_barplot_helper(tmp.counts, position, title, xlbl, ylbl)
    i <- i + 1
    min.n <- max.n + 1
  }
  return(out)
}

# Creates a single bar plot.
plot_barplot_helper <- function(wide.matrix, position, title, xlbl, ylbl, y.range = NULL,
                                labels = TRUE) {

  tmp <- reshape2::melt(t(wide.matrix)) # long format, correct order of categories
  names(tmp)[2] <- "Category"

  the.plot <- ggplot2::ggplot(tmp, ggplot2::aes(x = Var1, y = value, fill = Category))
  if (position == "stack") {
    the.plot <- the.plot + ggplot2::geom_col(pos = "stack")
  } else {
    the.plot <- the.plot + ggplot2::geom_col(pos = "fill")
  }
  if (!is.null(y.range)) {
    the.plot <- the.plot + ggplot2::coord_cartesian(ylim = y.range)
  }
  if (labels) {
    the.plot <- the.plot + ggplot2::labs(title = title, x = xlbl, y = ylbl) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))
  } else {
    the.plot <- the.plot +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
        legend.position = "none"
      )
  }
  return(the.plot)
}