% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_pca.R
\name{plot_PCA_ids}
\alias{plot_PCA_ids}
\title{Creates a PCA plot (by sample name).}
\usage{
plot_PCA_ids(PC.norm, PC.vars, components = 1:2)
}
\arguments{
\item{PC.norm}{Normalised PC matrix.}

\item{PC.vars}{Vector of variability explained by the PCs.}

\item{components}{The two components to include.}
}
\value{
A PCA plot.
}
\description{
Each sample is represented by their name.
}
\seealso{
\verb{\link{plot_PCA_ids_set}}

\verb{\link{plot_PCA_groups}}
}
