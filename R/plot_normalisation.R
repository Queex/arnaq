#' Creates a standard meanSdPlot
#' 
#' A wrapper for ggplot2's equivalent function.
#' @param assay Normalised data in a form `ggplot2::meanSdPlot` can interpret.
#' @returns The plot in question.
#' 
#' @export 
plot_meanSdPlot <- function(assay) {
  the.plot <- vsn::meanSdPlot(assay, plot = FALSE)$gg +
    ggplot2::labs(title = "Variance-stabilising normalisation")
  return(the.plot)
}