#' ARNAQ
#'
#' ARNAQ stands for *Analysts' RNA QC*. ARNAQ is a tool for exploring RNA sequencing datasets, with
#' an eye to differential expression or other analysis. It is designed to make reviewing the data
#' and isolating potential QC issues as simple and swift as possible.
#' 
#' It is not a replacement for the first line QC performed by sequencing platforms and their
#' associated tools; but instead to supplement them by considering experimental groups and allowing 
#' a more nuanced approach to potentially problematic samples than a simple pass/fail flag.
#'
#' Analysts who handle a lot of RNA projects will benefit from being able to perform these QC steps
#' more quickly, and will help avoid copy/paste scripts that are prone to errors where one line is
#' not updated for the current project. The report contains descriptions of what the plots show,
#' helping to present this work to collaborators who are less familiar with bioinformatics.
#'
#' Investigators who have some familiarity with bioinformatics will benefit from shorter, clearer R
#' scripts where they do not have to delve too deeply into the specifics of the packages used for
#' this QC analysis, helping them be confident that the analysis they perform is reliable.
#'
#' In both cases, the `.svg` export options the produces high quality print-ready plots that can be
#' easily incorporated into papers.
#'
#' *ARNAQ is named in reference to the board game 'Lost Ruins of Arnak'.*
#'
#' @docType package
#' @name arnaq
#' @seealso `\link{arnaq}`
NULL