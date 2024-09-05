#' Create a plot showing complexity
#' 
#' The complexity is measured by the proportion of the reads in the top n genes, as n varies, per
#' sample. Thinning is performed to reduce the file size for svg and html output).
#' 
#' To improve readability, this function can select the least complex samples to highlight, and
#' only show those in the legend. Non-highlighted samples will still have their curves drawn, just
#' in grey behind the highlighted samples.
#' 
#' In order to determine the 'worst' samples, a representative top n is used.
#' 
#' @param complexity.mat The matrix of calculated complexity, created automatically as part of the
#' main ARNAQ workflow.
#' @param highlight When `TRUE`, highlighting will be performed. When `FALSE`, all samples will be
#' drawn in colour, and will appear in the legend.
#' @param threshold The number of 'worst' samples to highlight, out of the total.
#' @param pivot What 'top n' to use when determining the worst samples.
#' @returns A plot showing proportion of reads assigned to the top n genes, for each sample, with
#' the least complex samples highlighted.
#' 
#' @export
plot_complexity_curves <- function(complexity.mat, highlight = TRUE, threshold = 8, pivot = 10) {
  if (highlight) {
    samp.levels <- levels(factor(complexity.mat$Sample))
    worst <- complexity.mat[complexity.mat$Count == pivot, ]
    worst <- worst[order(worst$Proportion, decreasing = TRUE), "Sample"]
    # Catch any samples with fewer than pivot genes
    sparse <- table(complexity.mat$Sample)
    sparse <- names(sparse)[sparse < pivot]
    worst <- c(sparse, as.character(worst))[seq_len(min(threshold, length(samp.levels)))]
    cols <- scales::hue_pal()(length(worst))
    names(cols) <- worst
  }

  # Thin the data to reduce draw time and filesize
  thinning.base <- sqrt(2)
  retain <- thinning.base^(0:floor(log(max(complexity.mat$Count), base = thinning.base)))
  retain <- round(retain)
  retain <- unique(c(
    retain,
    aggregate(complexity.mat$Count, by = list(complexity.mat$Sample), FUN = max)$x
  ))
  complexity.mat <- complexity.mat[complexity.mat$Count %in% retain, ]

  the.plot <- ggplot2::ggplot(
    complexity.mat,
    ggplot2::aes(x = Count, y = Proportion, colour = Sample, group = Sample)
  ) +
    ggplot2::geom_line() +
    ggplot2::scale_x_continuous(trans = "log2", breaks = c(1, 10, 100, 1000, 10000, 40000)) +
    ggplot2::scale_y_continuous(trans = "log2", breaks = c(5, 10, 20, 40, 60, 80, 100)) +
    ggplot2::labs(
      x = "Top n genes", y = "Percentage of reads for those genes",
      title = "Complexity curves"
    )

  if (highlight) {
    the.plot <- the.plot + ggplot2::scale_color_manual(values = cols, na.value = "#DFDFDF")
  }

  return(the.plot)
}

#' Create a plot showing detected genes
#' 
#' The plot shows how many genes are detected at certain count levels, to see how number of genes
#' drops off with higher counts.
#' 
#' To improve readability, this function can select the least complex samples to highlight, and
#' only show those in the legend. Non-highlighted samples will still have their curves drawn, just
#' in grey behind the highlighted samples.
#' 
#' In order to determine the 'worst' samples, a representative count is used.
#' 
#' @param masked.cpm.data The matrix of cpm data, masked to remove unwanted samples, created 
#' automatically as part of the
#' main ARNAQ workflow.
#' @param highlight When `TRUE`, highlighting will be performed. When `FALSE`, all samples will be
#' drawn in colour, and will appear in the legend.
#' @param threshold The number of 'worst' samples to highlight, out of the total.
#' @param pivot What count to use when determining the worst samples.
#' @returns A plot showing detected genes at different count levels, for each sample, with
#' the least complex samples highlighted.
#' @export
plot_detected_curves <- function(masked.cpm.data, highlight = TRUE, thresh = 8, pivot = 1) {

  cpm.values <- (1:100) / 10
  tmp <- reshape2::melt(apply(
    masked.cpm.data, 2,
    function(vec) {
      colSums(outer(vec, cpm.values, FUN = function(val1, val2) {
        val1 >= val2
      }))
    }
  ))
  colnames(tmp) <- c("CPM", "Sample", "Genes")
  tmp$CPM <- tmp$CPM / 10

  if (highlight) {
    samp.levels <- levels(tmp$Sample)
    worst <- tmp[tmp$CPM == pivot, ]
    worst <- worst[order(worst$Genes), "Sample"]
    # Catch any samples with fewer than pivot genes
    sparse <- table(tmp$Sample)
    sparse <- names(sparse)[sparse < pivot]
    worst <- c(sparse, as.character(worst))[seq_len(min(thresh, length(samp.levels)))]
    cols <- scales::hue_pal()(length(worst))
    names(cols) <- worst
  }

  the.plot <- ggplot2::ggplot(tmp,
                              ggplot2::aes(x = CPM, y = Genes, colour = Sample, group = Sample)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "CPM threshold", y = "Detected genes", title = "Detection curves")

  if (highlight) {
    the.plot <- the.plot + ggplot2::scale_color_manual(values = cols, na.value = "#DFDFDF")
  }

  return(the.plot)
}

# Creates the complexity matrix
make.complexity.matrix <- function(masked.data.mat, threshold = NULL,
                                   threshold.cpm = NULL) {
  cat("Calculating complexity matrix\n")

  if (is.null(threshold)) {
    threshold <- rep(0, ncol(masked.data.mat))
  }

  if (!is.null(threshold.cpm)) {
    # calculate on a per-sample basis
    threshold <- cbind(threshold, colSums(masked.data.mat) * threshold.cpm / 1000000)
    threshold <- apply(threshold, 1, max)
  }

  tmp.mat <- t(masked.data.mat)
  tmp.mat <- cbind(rownames(tmp.mat), threshold, tmp.mat)

  out <- apply(
    tmp.mat, 1,
    function(vec) {
      tmp <- make.complexity.curve(as.numeric(vec[c(-1, -2)]), as.numeric(vec[2]))
      if (length(tmp) > 0) {
        return(data.frame(Sample = vec[1], Count = seq_along(tmp), Proportion = tmp))
      } else {
        return(NULL)
      }
    }
  )
  out <- do.call("rbind", out)
  rownames(out) <- NULL
  return(out)
}

# Helper function to create the curve for a single sample for make.complexity.matrix
make.complexity.curve <- function(masked.data.vec, threshold = 0) {
  total <- sum(masked.data.vec) # Get divisor
  tmp <- masked.data.vec[masked.data.vec > threshold] # Get only present genes
  # Put genes in descending order of read depth and build cumulative function
  tmp <- cumsum(sort(tmp, decreasing = TRUE)) / total
  return(tmp * 100) # turn into percentage
}

# Creates the cpm matrix
make.cpm.matrix <- function(masked.data.mat) {
  cat("Calculating detection matrix\n")
  return(t(t(masked.data.mat) * 1000000 / colSums(masked.data.mat)))
}