# Helper function to recursively export the `all_plots` list.
export_all_svg.helper <- function(name.fragment, object) {
  if (ggplot2::is.ggplot(object)) {
    export_svg(object, paste0(name.fragment, ".svg"))
  } else { # need to recurse
    if (is.null(names(object))) { # by number
      for (i in seq_along(object)) {
        export_all_svg.helper(paste(name.fragment, i, sep = "_"), object[[i]])
      }
    } else { # by name
      for (n in names(object)) {
        export_all_svg.helper(paste(name.fragment, n, sep = "_"), object[[n]])
      }
    }
  }
}

#' Exports all the plots created by arnaq in svg format.
#' 
#' This function is called automatically when `\link{arnaq}` is called with `svg=TRUE`. However,
#' it can be called after adjusting any of the generated plots to ensure the svg export is of the
#' correct versions of them.
#' @param svg.directory The directory into which to make the export. If it does not exist, it will 
#' be created, if possible. The default is `QC/svg`.
#' @param qc.name The prefix to use for all the file names. By default, this is the project name -
#' underscore - model name as set in `resources.yml` and the `model.name` argument to 
#' `\link{arnaq}`.
#' @seealso `\link{export_svg}`
#' @seealso `\link{arnaq}`
#' 
#' @export
export_all_svg <- function(svg.directory = "QC/svg", qc.name = qc.name) {
  cat("\nSaving svg plots\n")
  ensure.directory(svg.directory)
  export_all_svg.helper(paste0(svg.directory, "/", qc.name), all_plots)
  cat("Saving svg plots done\n\n")
}

#' Exports a single plot as svg.
#' 
#' A convenience function for exporting svg plots, in a consistent format with bulk exports from
#' `\link{export_all_svg}`.
#' @param the.plot The plot to export.
#' @param filename The filename to export to.
#' @seealso `\link{export_all_svg}`
#' 
#' @export
export_svg <- function(the.plot, filename) {
  cat(paste("Saving svg to", filename, "\n"))
  ggplot2::ggsave(file = filename, plot = the.plot, width = 6, height = 6)
}