#' Create ARNAQ support files
#'
#' Creates copies of important ARNAQ files in the working directory.
#'
#' `arnaq_template_1.0.rmd` is a skeleton rmarkdown object used to generate the report. Besides
#' using the default version, you can edit to to include branding, citations for the specific
#' platform and tools used, etc..
#'
#' `biotype_conversion_table.txt` is used to map the full set of Ensembl biotypes to a more
#' manageable set for plotting. The default version should be fine for most applications, but if you
#' want finer control over how the categories are collapsed you can create your own based on the
#' default.
#'
#' `resources.yml` is the file used to give project-wide metadata, and an example file is provided.
#'
#' @param files A vector naming which support files you want copies of. The default consists of
#' all available files. Any unrecognised entries will be ignored.
#' @param ... Further arguments to be passed to file.copy (such as `overwrite=TRUE`).
#'
#' @export
arnaq.prepare <- function(files=c("resources", "template", "biotypes"), ...) {
  if ("resources" %in% files) {
    arnaq.prepare.helper("resources.yml", ...)
  }
  if ("template" %in% files) {
    arnaq.prepare.helper("arnaq_template_1.0.rmd", ...)
  }
  if ("biotypes" %in% files) {
    arnaq.prepare.helper("biotype_conversion_table.txt", ...)
  }
  invisible()
}

arnaq.prepare.helper <- function(file, ...) {
  cat(paste("Creating", file, "in working directory.\n"))
  from.file <- system.file("extdata", file, package="arnaq")
  file.copy(from.file, "./", ...)
}