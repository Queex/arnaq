#' Creates a plot of observed and expected ERCC spike-in concentrations
#' 
#' This is intended where a single ERCC mix has been used across multiple samples. The observed read
#' counts are restricted to the ERCC spike-ins and concentration estimates constructed from these
#' values.
#' 
#' The `samples.txt` file (and hence `sample.metadata`) should contain ERCC mix data when this 
#' function is
#' used; a column named `ERCC` containing the values 1 and 2 to specify which mix.
#' 
#' @param read_matrix The raw count matrix, samples as columns.
#' @param sample.metadata The sample metadata.
#' @param mix Either 1 or 2, to indicate which ERCC mix was used.
#' @param ERCC.ref.table The table of ERCC spike-in data.
#' @returns A plot showing observed ERCC spike-in proportions against expected.
#' 
#' @export
plot_ERCC_observed_scatter <- function(read_matrix, sample.metadata, mix,
                                       ERCC.ref.table = ERCC.table) {

  tmp.ERCC.table <- ERCC.ref.table[match(rownames(read_matrix), ERCC.ref.table[, 2]), ]
  expected <- tmp.ERCC.table[, 3 + mix]

  dat <- read_matrix[, sample.metadata$ERCC == paste0("ERCC", mix), drop = FALSE]
  dat <- reshape2::melt(t(t(dat) / colSums(dat)))
  names(dat) <- c("ERCC", "Sample", "Observed")
  dat$Expected <- expected / sum(expected)
  dat <- dat[dat$Observed > 0, ] # Remove 0 observed

  dat.corr <- round(stats::cor(dat$Expected, dat$Observed), 3)

  the.plot <- ggplot2::ggplot(dat, ggplot2::aes(x = Expected, y = Observed, colour = Sample)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "glm", colour = "red") +
    ggplot2::geom_abline(intercept = 0, slope = 1, colour = "green") +
    ggplot2::scale_x_continuous(trans = "log2",
                                breaks = c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1)) +
    ggplot2::scale_y_continuous(trans = "log2",
                                breaks = c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1)) +
    ggplot2::ggtitle(paste("Observed and expected Spike-ins for mix", mix)) +
    ggplot2::labs(y = "Observed CPM", x = "Expected CPM") +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(caption = paste("Correlation:", dat.corr))

  return(the.plot)
}

#' Creates a plot of errors in fold-change estimates for ERCC spike-ins
#' 
#' The fold-change is calculated between two specific samples, and all such pairs are summarised
#' in this plot, showing the estimate errors against the true values given by the eRCC spike-in 
#' definitions.
#' 
#' @param mix1_reads A dataframe of normalised reads for the mix 1 samples.
#' @param mix2_reads A dataframe of normalised reads for the mix 2 samples.
#' @param ERCC_names A vector of names of the ERCC spike-ins
#' @param title The title to give the plot
#' @param ERCC.ref.table The table of ERCC spike-in data.
#' @returns A plot showing the estimates in fold-change
#' @seealso `\link{plot_all.ERCC.pairs}`
#' 
#' @export
plot_ERCC_errors <- function(mix1_reads, mix2_reads, ERCC_names, title,
                             ERCC.ref.table = ERCC.table) {

  logratio <- as.vector(log2(mix1_reads / mix2_reads))
  # Need to sort ERCC table
  tmp.ERCC.table <- ERCC.ref.table[match(ERCC_names, ERCC.ref.table[, 2]), ]
  net_concentration <- sqrt(tmp.ERCC.table[, 4] * tmp.ERCC.table[, 5])
  # geomean for sensible consensus value
  expected <- factor(tmp.ERCC.table[, 7])

  # Make logratio the absolute error
  logratio <- logratio - as.numeric(as.character(expected))

  palette <- scales::hue_pal()(4)
  names(palette) <- c("0", "2", "-1", "-0.58")

  dat <- data.frame(
    Observed = logratio, Expected = expected,
    Log10Concentration = log10(net_concentration)
  )
  dat$Detected <- is.finite(dat$Observed)
  dat$Observed[!dat$Detected] <- 0
  dat$Detected <- factor(as.character(dat$Detected), levels = c("TRUE", "FALSE"))

  the.plot <- ggplot2::ggplot(dat,
    ggplot2::aes(x = Log10Concentration, y = Observed, colour = Expected, shape = Detected)
  ) +
    ggplot2::scale_shape_manual(values = c(19, 4)) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::ggtitle(title) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::labs(x = "Log10 Concentration", y = "Error in f/c estimate, log2 scale") +
    ggplot2::geom_point()

  return(the.plot)
}

#' Creates a scatterplot of ERCC spike-in fold-change ratios against the true values
#' 
#' Takes two samples, one from each mix, and calculates the fold change between them for ERCC 
#' spike-ins, the creates a scatter plot against the true, known fold change values.
#' @param read_matrix The full normalised count matrix.
#' @param mix1_column The name of the sample with ERCC mix 1.
#' @param mix2_column The name of the sample with ERCC mix 2.
#' @param ERCC.ref.table The table of ERCC spike-in data.
#' @returns A scatter plot of the fold-change against exected values.
#' @seealso `\link{plot_all.ERCC.pairs}`
#' 
#' @export
plot_ERCC_fc_scatter <- function(read_matrix, mix1_column, mix2_column,
                                 ERCC.ref.table = ERCC.table) {

  logratio <- log2(read_matrix[, mix1_column] / read_matrix[, mix2_column])

  # Need to sort ERCC table
  tmp.ERCC.table <- ERCC.ref.table[match(rownames(read_matrix), ERCC.ref.table[, 2]), ]
  net_concentration <- sqrt(tmp.ERCC.table[, 4] * tmp.ERCC.table[, 5])
  # geomean for sensible consensus value

  expected <- tmp.ERCC.table[, 7]

  num.nan <- sum(is.nan(logratio))
  num.inf <- sum(!is.finite(logratio))

  # find correlation
  not.nan.mask <- !is.nan(logratio)
  all.quantity.mask <- is.finite(logratio)
  high.quantity.mask <- all.quantity.mask & (net_concentration > stats::median(net_concentration))

  high.quantity.cor <- round(stats::cor(logratio[high.quantity.mask],
                                        expected[high.quantity.mask]), 3)
  all.quantity.cor <- round(stats::cor(logratio[all.quantity.mask], expected[all.quantity.mask]), 3)

  # make table
  dat <- data.frame(Observed = logratio, Expected = expected, Concentration = net_concentration)
  dat <- dat[all.quantity.mask, ]

  theplot <- ggplot2::ggplot(dat, 
                             ggplot2::aes(x = Expected, y = Observed, colour = Concentration)) +
    ggplot2::geom_jitter(shape = 1, height = 0, width = 0.05) +
    ggplot2::scale_colour_gradient(low = "#F0F0F0", high = "#000000", trans = "log") +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "red") +
    ggplot2::ggtitle(paste0("ERCC fold-changes\n", mix1_column, " to ", mix2_column)) +
    ggplot2::labs(caption = paste0(
      "Below detection: ",
      num.nan + num.inf,
      "  Correlation: high-quantity=",
      high.quantity.cor, " all=", all.quantity.cor
    ))

  return(theplot)
}

#' Creates plots showing ERCC spike-in comparisons
#' 
#' The behaviour of this function is controlled by `ERCC.combined`. When `TRUE`, all the pairs of
#' samples are combined into a single plot showing the fold change estimate errors. Otherwise, it
#' creates a list of plots, one for each pair, showing the scatterplot of observed against
#' expected fold-changes.
#' 
#' @param ERCC.data.cpm The ERCC data, normalised to counts per million format.
#' @param ERCC.pairs A list of length 2 vectors of sample names, each pair representing one 
#' comparison.
#' @param ERCC.combined If `TRUE`, consolidate the pairs into a singple plot, otherwise produce a 
#' series of separate plots.
#' @returns A list containing one or more plots.
#' @seealso `\link{plot_ERCC_observed_scatter}`
#' @seealso `\link{plot_ERCC_errors}`
#' @seealso `\link{plot_ERCC_fc_scatter}`
#' 
#' @export
plot_all.ERCC.pairs <- function(ERCC.data.cpm, ERCC.pairs, ERCC.combined) {
  out <- list()
  if (ERCC.combined) {
    ERCC_names <- rownames(ERCC.data.cpm)
    if (length(ERCC.pairs) > 1) { # Make combined plot if >1 ERCC comparison
      tmp1 <- reshape2::melt(ERCC.data.cpm[, unlist(lapply(ERCC.pairs, function(x) {
        return(x[1])
      }))])
      tmp2 <- reshape2::melt(ERCC.data.cpm[, unlist(lapply(ERCC.pairs, function(x) {
        return(x[2])
      }))])
      tmp <- plot_ERCC_errors(
        tmp1[, 3], tmp2[, 3],
        rep(ERCC_names, length(ERCC.pairs)),
        paste0("ERCC fold-changes\nAll pairs")
      )
      out[["All"]] <- tmp
    }
    for (pair in ERCC.pairs) {
      tmp <- plot_ERCC_errors(
        ERCC.data.cpm[, pair[1], drop = FALSE],
        ERCC.data.cpm[, pair[2], drop = FALSE], ERCC_names,
        paste0("ERCC fold-changes\n", pair[1], " to ", pair[2])
      )
      out[[paste(pair[1], "v", pair[2])]] <- tmp
    }
  } else {
    for (pair in ERCC.pairs) {
      tmp <- plot_ERCC_fc_scatter(ERCC.data.cpm, pair[1], pair[2])
      out[[paste(pair[1], "v", pair[2])]] <- tmp
    }
  }
  return(out)
}

# Helper function to create log2 fold-change values
make.fc.dat <- function(ERCC.dat, pair) {
  return(log2(ERCC.dat[, pair[1]] / ERCC.dat[, pair[2]]))
}