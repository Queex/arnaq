######################################################
# Analyst's RNA QC
# v 0.0.1
######################################################
# Main entry point:
# arnaq                  arnaq()
######################################################
# TODO
# Add heatmap/dendro/PCA? based on most variable
# Utility functions for showing what plots have been made
# Make log prettier, include timestamps
# Improve plots (especially dendro)
# Infer plot flags from plot object
# Check that Picard's duplication metrics file works as-is
# Adds checks for sample order to metrics files
# Change violin plots to log scale
# Automatically remove Groups with a single level if treat.groups is not specified
######################################################

#' Create an ARNAQ Report
#'
#' This function creates a report based on the RNA count data, genome resources and metadata
#' resources, in html format. It uses external files to define specific file locations and
#' parameters rather than incorporating them all into the function arguments.
#'
#' `resources.yml` contains information that is identical for all samples, or refers to the
#' project
#' as a whole.
#'
#' `samples.txt` contains sample-centric information, such as membership of an experimental
#' group.
#'
#' Both of these are assumed to remain unchanged through different ARNAQ runs. Arguments to this
#' function are those options that may change run to run. More details are in ????.
#'
#' A directory will be created in the current working directory called
#' `QC` to house the output files, over-writing previous versions, if any.
#'
#' As a side effect of running this function, a number of R objects are created in the session.
#' The intention is that these are available for downstream analysis, but can also be adjusted
#' (such as fine-tuning plots). The R console output will also appear in the file `arnaq.log`
#' in your working directory.
#'
#' ARNAQ will try to optimise by looking in the session for some data before attempting to read
#' it from the file locations you give it. If the underlying files have changed, use
#' `\link{arnaq_clear}` to remove ARNAQ's objects before rerunning.
#'
#' @param resources.file The file containing project-centric data, including the location of counts
#' tables and metric tables. The default looks for the file `resources.yml` in your working
#' directory, and is useful when you
#' have one project per directory.
#' @param samples.txt The file containin sample-centric data. The default looks for the file
#' `samples.txt` in your working directory, and is useful when you have one project per
#' directory.
#' @param model.name An optional model name that allows you to have multiple ARNAQ reports in
#' parallel for the same project. The model name is appended to the project name in ARNAQ outputs.
#' Note that running ARNAQ will replace its generated data objects in the current R session.
#' @param sample.mask A vector of True/False values, with length equal to the number of samples,
#' indicating which samples to include in the report. This is used to run reports excluding
#' outliers. `\link{make_outlier_mask}` is a convenient function for generating one. Default
#' behaviour is to include all samples.
#' @param treat.groups A vector of character objects that specifices which of the group columns in
#' the `samples.txt` file are to be considered treatment groups for `DESeq2`'s
#' normalisation. If the group information is not of full rank, this should be a constrained set of
#' these groups that is of full rank. Default behaviour is to use all groups.
#' @param pca.groups A vector of character objects that specify which groups in `samples.txt`
#' to create PCA plots for. Default behaviour
#' is to only use those groups in `treat.groups`.
#' @param ERCC A logical value indicating if ERCC spike-ins are included in the gene set.
#' @param normalise A character value specifying which normalisation to perform on the data.
#' `VSN` is the default, but this does not run in reaosnable time for large (>=30) numbers
#' of samples. `linear` is the alternative. Any other value will skip normalisation.
#' @param svg.export If set to `TRUE`, create an additional directory in the output called
#' `svg` with svg versions of every plot created in the report, suitable for inclusion in
#' papers.
#' @param max.sample.per.page In projects with a large number of samples, some plots will be
#' spread across multiple 'pages' in the output report, to keep them readable. Alter the default
#' value to fine-tune the point at which this behaviour will occur.
#' @param pca.depth How many components to use when preparing PCA plots. Components beyond 2 will
#' be plotted against the first component, rather than showing every pairwise combination up to the
#' limit. Default behaviour is to use all components that account for at least 10% of the variation,
#' with a minimum of 2.
#' @param ERCC.pairs A list of pairs of samples to produce ERCC-specific fold-change plots for.
#' Each entry in the list should be a length 2 vector of character objects, naming those samples.
#' Should be an empty list if ERCC spike-ins are not present in the data.
#' @param ERCC.combined When `TRUE`, the ERCC.pairs plots above will use concentration as the
#' x-axis, rather than the expected fold-change.
#' @param scatter.pairs A list of pairs of samples to produce read counts scatter plots for, to
#' check consistency between samples that are technical replicates.
#' @param gene.mask.name Which set of genes to filter to in preparing the report. This should be a
#' name that will be present in the `gene.masks` list. See the `objects created` vignette for a list
#' of options. The default of `Genes` will use all genes in the reference `.gtf` with a non-zero
#' count for at least one sample, that are not ERCC spike-ins.
#' @seealso \code{\link{make_outlier_mask}}
#' @seealso `\link{arnaq_clear}`
#' @examples
#' arnaq()
#' arnaq("some/other/location/resources.yml", "some/other/location/samples.txt")
#' arnaq(model.name = "outliers removed", sample.mask = some_mask, normalise = "linear")
#'
#' @export
arnaq <- function(resources.file = "resources.yml",
                  sample.file = "samples.txt", model.name = NULL, sample.mask = NULL,
                  treat.groups = NULL, pca.groups = NULL, ERCC = FALSE, normalise = "VSN",
                  svg.export = FALSE, max.samples.per.page = 40, pca.depth = NULL,
                  ERCC.pairs = list(), ERCC.combined = TRUE, scatter.pairs = list(),
                  gene.mask.name = "Genes") {

  clear.sinks()
  sink("arnaq.log", split = TRUE)

  arnaq.version <- "0.1"
  template.version <- "1.0"

  cat("Starting ARNAQ\n")
  cat(paste("Script version:", arnaq.version, "\n\n"))
  resources <- read.resources.file(resources.file)
  required.resources <- c(
    "project_id", "count_table", "report_template",
    "resource_dir", "species"
  )
  optional.resources <- c(
    "summary_table", "duplication_table", "metrics_table",
    "ercc_concentrations", "biotype_conversion", "genome_reference"
  )
  if (!all(names(resources) %in% c(required.resources, optional.resources))) {
    warning("Unexpected line in resource file!")
    warning(names(resources)[! names(resources) %in% c(required.resources, optional.resources)])
  }
  if (!all(c(required.resources, optional.resources) %in% names(resources))) {
    stop("Not all expected lines present in resources file!")
  }
  for (n in required.resources) {
    if (resources[[n]] == "????" || resources[[n]] == "None" || resources[[n]] == "none") {
      stop(paste("Missing value for required resource in resource file:", n))
    }
  }

  arnaq.run <- list(version=arnaq.version, date=Sys.Date(), resources=resources)

  arnaq.run$project.id <- resources[["project_id"]]
  count.table <- resources[["count_table"]]
  assignment.table <- resources[["summary_table"]]
  duprate.table <- resources[["duplication_table"]]
  metrics.table <- resources[["metrics_table"]]
  ercc.table <- resources[["ercc_concentrations"]]
  resource.dir <- resources[["resource_dir"]]
  species <- resources[["species"]]
  biotype.conversion <- resources[["biotype_conversion"]]
  if (biotype.conversion == "INTERNAL") {
    biotype.conversion <- system.file("extdata", "biotype_conversion_table.txt", package="arnaq")
  }
  arnaq.report.template <- resources[["report_template"]]
  if (arnaq.report.template == "INTERNAL") {
    arnaq.report.template <- system.file("extdata",
                                         paste0("arnaq_template_", template.version, ".rmd"),
                                         package="arnaq")
  }
  arnaq.report.template <<- arnaq.report.template
  arnaq.run$genome.file <- resources[["genome_reference"]]
  expected.template.version <- "0.1"
  arnaq.run$out.directory <- "QC/"

  check.QC.template(arnaq.report.template, expected.template.version)

  cat("\n")

  qc.name <<- ifelse(is.null(model.name), arnaq.run$project.id,
    paste(arnaq.run$project.id, model.name, sep = "_")
  )
  cat(paste("This QC run is named:", qc.name, "\n"))

  # Set up output directory
  ensure.directory(arnaq.run$out.directory)

  # Load data
  count.data <<- read.count.data(c("count.data", count.table))

  # Read sample metadata
  sample.metadata <<- read.samples.metadata(c("sample.metadata", sample.file), count.data)

  # Rename count.data columns
  colnames(count.data) <<- sample.metadata$Display

  # Load gtf if available
  if (arnaq.run$genome.file == "None" || arnaq.run$genome.file == "none") {
    species.gtf <<- NULL
  } else {
    species.gtf <<- read.biotypes(arnaq.run$genome.file, resource.dir)
  }

  # Make gene masks
  gene.masks <<- make.gene.masks(count.data, species.gtf, ERCC = ERCC)

  featureCount.metrics <- assignment.table != "None"
  if (featureCount.metrics) {
    arnaq.run$read.summary <- read.read.summary(c("arnaq.run$read.summary", assignment.table),
                                                count.data, gene.masks, sample.metadata)
    colnames(arnaq.run$read.summary) <- sample.metadata$Display
  }

  # Read extra metrics
  picard.metrics <- (metrics.table != "None") && (metrics.table != "none")
  if (picard.metrics) {
    arnaq.run$count.metrics <- read.picard.metrics(metrics.table, duprate.table,
                                                   sample.metadata$Display)
  }

  # Sample mask
  if (is.null(sample.mask)) {
    sample.mask <- TRUE
    outliers.removed <- FALSE
  } else {
    cat("Applying outlier mask.\n")
    outliers.removed <- TRUE
  }

  assign("sample.mask", sample.mask, 1)

  # Check that samples named in scatter and ERCC exist in data
  tmp <- unique(c(unlist(ERCC.pairs), unlist(scatter.pairs)))
  tmp2 <- !(tmp %in% sample.metadata$Display[sample.mask])
  if (sum(tmp2) > 0) {
    stop(paste(
      "Some named samples are not present in the data:",
      paste(tmp[tmp2], collapse = " ")
    ))
  }

  # Adjust max samples per page to be prettier
  n <- sum(sample.mask)
  if (n > max.samples.per.page) {
    pages <- ceiling(n / max.samples.per.page)
    max.samples.per.page <- ceiling(n / pages)
  }

  # Make temp data
  tmp.sample.metadata <- sample.metadata[sample.mask, , drop=FALSE]
  tmp.count.data <- count.data[gene.masks[[gene.mask.name]], sample.mask, drop=FALSE]
  if (featureCount.metrics) {
    tmp.read.summary <- arnaq.run$read.summary[, sample.mask, drop=FALSE]
  } else {
    tmp.read.summary <- list()
  }
  if (picard.metrics) {
    tmp.count.metrics <- list(
      ByBase = arnaq.run$count.metrics$ByBase[sample.mask, , drop=FALSE],
      ByOther = arnaq.run$count.metrics$ByOther[sample.mask, , drop=FALSE]
    )
  } else {
    tmp.count.metrics <- list()
  }

  # Export filtered sample metadata
  save.counts(tmp.sample.metadata,
              paste0(arnaq.run$out.directory, "/", qc.name, "_sample_metadata.txt"))

  if (outliers.removed) {
    cat(paste("Remaining samples:", paste(tmp.sample.metadata$Display, collapse = " "), "\n"))
  }

  # Save counts
  save.counts(tmp.count.data,
              paste0(arnaq.run$out.directory, "/", qc.name, "_gene_count_table.txt"))
  if (ERCC) {
    read.ERCC.table(ercc.table)
    arnaq.run$ERCC.data <- count.data[gene.masks$ERCC, sample.mask, drop=FALSE]
    if (any(colSums(arnaq.run$ERCC.data) == 0)) {
      warning("At least one sample has no ERCC counts.")
    }
    save.counts(arnaq.run$ERCC.data,
                paste0(arnaq.run$out.directory, "/", qc.name, "_ERCC_count_table.txt"))
  }

  # Save detected genes
  tmp.detected <- cbind(Genes = colSums(tmp.count.data > 0))
  if (ERCC) {
    tmp.detected <- cbind(tmp.detected, ERCC = colSums(arnaq.run$ERCC.data > 0))
  }
  save.counts(tmp.detected,
              paste0(arnaq.run$out.directory, "/", qc.name, "_detected_genes_table.txt"))

  # Define treatment groups
  if (is.null(treat.groups)) {
    treat.cols <- 4:ncol(tmp.sample.metadata)
  } else {
    treat.cols <-
      (seq_len(ncol(tmp.sample.metadata)))[colnames(tmp.sample.metadata) %in% treat.groups]
  }
  primary.treat.group <- treat.cols[1]
  cat(
    "Considering these columns as Treatments:",
    paste(colnames(tmp.sample.metadata[treat.cols]), collapse = " "), "\n"
  )
  cat("Primary Treatment group:", colnames(tmp.sample.metadata)[primary.treat.group], "\n")

  # Define columns for PCA
  if (is.null(pca.groups)) {
    pca.cols <- treat.cols
  } else {
    pca.cols <- (seq_len(ncol(tmp.sample.metadata)))[colnames(tmp.sample.metadata) %in% pca.groups]
  }
  cat(
    "These columns will be displayed in the PCA plots:",
    paste(colnames(tmp.sample.metadata[pca.cols]), collapse = " "), "\n"
  )

  cat("\n")

  # Calculate complexity
  arnaq.run$complexity.matrix <- make.complexity.matrix(tmp.count.data, threshold.cpm = 1)

  # Calculate cpm
  arnaq.run$cpm.matrix <- make.cpm.matrix(tmp.count.data)
  save.counts(arnaq.run$cpm.matrix, paste0(arnaq.run$out.directory, "/", qc.name, "_cpm_table.txt"))

  # Use DESeq2
  if (normalise == "VSN" || normalise == "vsn") {
    cat("Applying variance stabilising normalisation\n")
    arnaq.run$model.formula <- formula(paste0(
      "~",
      paste(colnames(tmp.sample.metadata)[treat.cols], collapse = " + ")
    ))
    tmp <- apply.DESeq2(tmp.count.data, tmp.sample.metadata, arnaq.run$model.formula)
    cds.deseq2 <- tmp[[1]]
    norm.data <<- tmp[[2]] # Blinded normalisation to use for any downstream analysis.

    # Save normalised counts
    save.counts(norm.data,
                paste0(arnaq.run$out.directory, "/", qc.name, "_normalised_count_table.txt"))
  } else if (normalise == "linear" || normalise == "Linear") {
    cat("Applying linear normalisation\n")
    norm.data <<- t(t(tmp.count.data) / colSums(tmp.count.data))
  } else {
    cat("Leaving data unnormalised\n")
    norm.data <<- tmp.count.data
  }

  # Calculate PCA
  arnaq.run$norm.PCA <- calc.pca(norm.data)
  if (is.null(pca.depth)) { # calculate how many based on 10% threshold
    pca.depth <- which(arnaq.run$norm.PCA$Vars < 10)[1] - 1
    pca.depth <- max(2, pca.depth)
  }

  # Simple normalisation for ERCC data
  if (ERCC) {
    arnaq.run$ERCC.data.cpm <- apply(arnaq.run$ERCC.data, 2, function(vec) {
      vec * 1000000 / sum(vec)
    })
  }

  # Create biotype data, if available
  if (!is.null(species.gtf)) {
    biotype.conversion.table <- read.biotype.conversion.table(biotype.conversion)
    arnaq.run$biotype.conversion.table <- biotype.conversion.table
    arnaq.run$biotype.counts <- make.biotype.data(tmp.count.data,
                                                  species.gtf, biotype.conversion.table)
    do.biotype <- length(arnaq.run$biotype.counts) > 0
  } else {
    do.biotype <- FALSE
  }
  if (do.biotype) {
    save.counts(arnaq.run$biotype.counts[[1]],
                paste0(arnaq.run$out.directory, "/", qc.name, "_biotypes_table.txt"))
    save.counts(arnaq.run$biotype.counts[[3]],
                paste0(arnaq.run$out.directory, "/", qc.name, "_biotypes_original_table.txt"))
  }

  # Assign arnaq.run to .GlobalEnv to guarantee it is not lost
  assign("arnaq.run", arnaq.run, 1)

  # Plot
  arnaq_report(arnaq.run,
    arnaq.report.template, qc.name, arnaq.run$genome.file, norm.data, tmp.count.data,
    tmp.sample.metadata, primary.treat.group, tmp.read.summary, tmp.count.metrics,
    arnaq.run$biotype.counts, arnaq.run$complexity.matrix, arnaq.run$cpm.matrix,
    arnaq.run$out.directory, svg.export, ERCC, featureCount.metrics,
    picard.metrics, normalise, do.biotype, max.samples.per.page, pca.depth, pca.cols,
    ERCC.pairs, ERCC.combined, scatter.pairs
  )

  cat("ARNAQ finished!\n")

  # Close log
  sink(NULL)
}

#' Create a mask based on sample names
#'
#' Returns a logical vector for the samples in the project that do not match any name in the
#' supplied vector. This is intended for use with the `sample.mask` parameter of the
#' `\link{arnaq}` function. The names are compared to the Display column of `samples.txt`.
#'
#' @param outlier.names A vector of character objects with the names of the samples to
#' \strong{exclude}.
#' @returns A logical vector, `TRUE` for samples to include.
#' @seealso `\link{arnaq}`
#'
#' @export
make_outlier_mask <- function(outlier.names) {
  # Check that all given names are in sample.metadata
  if (sum(!(outlier.names %in% sample.metadata$Display)) > 0) {
    warning("You gave a sample name not present in the data")
  }

  # return mask
  !(sample.metadata$Display %in% outlier.names)
}

#' Removes ARNAQ objects
#'
#' Removes objects ARNAQ creates from the current session. This is useful to force ARNAQ to load
#' files from disk if they have been changed, or to 'clean up' the session for downstream analysis.
#' Any warnings (such as due to the objects not existing) are suppressed.
#'
#' It is good practice to put `arnaq_clear()` at the top of your script so if you rerun from scratch
#' you can be sure the data is created afresh.
#'
#' @param all If the default of `TRUE`, all of ARNAQ's objects will be removed. Otherwise,
#' `count.data`, `gene.masks`, `sample.mask` and `sample.metadata` will be kept, as the objects
#' needed for downstream analysis.
#'
#' @export
arnaq_clear <- function(remove.all=TRUE) {
  always.remove <- c("all_plots", "arnaq.plot.flags", "arnaq.report.template", "arnaq.run",
                     "qc.name", "norm.data", "species.gtf")
  optional.remove <- c("count.data", "gene.masks", "sample.mask", "sample.metadata")
  if (remove.all) {
    always.remove <- c(always.remove, optional.remove)
  }
  suppressWarnings(rm(list=always.remove, pos = 1))
}