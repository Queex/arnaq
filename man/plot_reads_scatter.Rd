% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_scatter.R
\name{plot_reads_scatter}
\alias{plot_reads_scatter}
\title{Creates a scatter plot of reads for two samples.}
\usage{
plot_reads_scatter(
  norm.dat,
  dat.for.filter = norm.dat,
  col1name,
  col2name,
  title = NULL,
  log = FALSE,
  filter.threshold = 3
)
}
\arguments{
\item{norm.dat}{The dataframe containing the reads to base the scatterplot on.}

\item{dat.for.filter}{A similar dataframe containing the data to use for filtering. By default
the same as \code{norm.dat}.}

\item{col1name}{The name of the coolumn containing data for sample 1.}

\item{col2name}{The name of the coolumn containing data for sample 2.}

\item{title}{The title to use for the plot.}

\item{log}{If \code{TRUE}, use log2 scales for the axes.}

\item{filter.threshold}{The minimum value to filter on in \code{dat.for.filter}. The default of 3
prevents numerical issues when preparing the plot.}
}
\value{
a scatterplot.
}
\description{
The scatter plot compares the reads for genes across two samples, preferably using normalised
data. A filter is applied for minimum read depth for a gene, optionally using a different
dat matrix (e.g. unnormalised).
}
\seealso{
\verb{\link{plot_all_scatter_pairs}}
}
