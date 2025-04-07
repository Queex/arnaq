# Creates all applicable plots and generates final report
arnaq_report <- function(arnaq.run = arnaq.run,
                         rmarkdown_template, qc.name, genome.file, norm.data, count.data,
                         sample.metadata, primary.treat.group, read.summary, tmp.count.metrics,
                         biotype.counts, complexity.mat,
                         cpm.mat, plot_directory, svg.export, ERCC = TRUE,
                         featureCount.metrics = TRUE, picard.metrics = TRUE,
                         normalised = "VSN", do.biotype, max.samples.per.page = 40,
                         pca.depth = 2, pca.cols = 4,
                         ERCC.pairs = list(), ERCC.combined = TRUE, scatter.pairs = list()) {

  # Make sure tinytex is installed - this is not done when installing the package but requires
  # a call from within R.
  # If there are tex-related problems, running tinytex::uninstall() before trying this
  # script again might fix it.
  if (!tinytex::is_tinytex()) {
    tinytex::install_tinytex(force = TRUE) # Suppress 'existing installation' warning
  }

  n.samples <- ncol(count.data)

  # Boolean flags for different sections
  # Stop showing fully-labelled complexity plot when there are too many samples
  complexity_plot <- (n.samples <= 40)
  # Show 'highlights' complexity plot once there are more than 12 samples
  complexity_plot_2 <- (n.samples > 12)
  ERCC_mix <- (ERCC && "ERCC" %in% names(sample.metadata))
  ERCC_fc <- (ERCC && length(ERCC.pairs) > 0)
  ERCC_show <- (ERCC_mix || ERCC_fc)
  scatterplots <- (length(scatter.pairs) > 0)

  all_plots <- list()

  # barplots
  if (featureCount.metrics) {
    cat("Making barplots\n")
    all_plots$barplots_ps <- plot_counts_barplot(read.summary, , max.samples.per.page)
    all_plots$barplots2_ps <- plot_counts_barplot(read.summary, , max.samples.per.page, prop = TRUE)
  }

  # violin plots of extra metrics
  if (picard.metrics) {
    tmp <- plot_all_metrics_violin(tmp.count.metrics, primary.treat.group)
    all_plots$violin1 <- tmp[[1]]
    all_plots$violin2 <- tmp[[2]]
  }

  # Transformation plots
  if (normalised == "VSN") {
    cat("Making normalisation plot\n")
    all_plots$norm_plot <- plot_meanSdPlot(as.matrix(norm.data))
    # The difference between this and norm.data is that this method is unblinded
    # with respect to the experimental design. Better for analysing the quality of the data but
    # should not be used for downstream analysis as it would introduce bias.
  }

  # 4 hierarchical clustering and heatmap plots
  cat("Making hc plots\n")
  all_plots$hcluster_ps <- plot_all_clustering(norm.data)

  # PCA
  cat("Making PCA plots\n")
  all_plots$pca_ids_ps <- plot_PCA_ids_set(arnaq.run$norm.PCA$PC, arnaq.run$norm.PCA$Vars,
                                           depth = pca.depth)
  all_plots$pca_group_pss <- plot_all_PCA_groups(arnaq.run$norm.PCA$PC, arnaq.run$norm.PCA$Vars,
                                                 sample.metadata, pca.cols, depth = pca.depth)

  # Complexity
  cat("Making complexity plot(s)\n")
  if (complexity_plot) {
    all_plots$complexity_plot <- plot_complexity_curves(complexity.mat, FALSE)
  }
  if (complexity_plot_2) {
    all_plots$complexity_plot_2 <- plot_complexity_curves(complexity.mat, TRUE)
  }

  # Detection
  cat("Making detected genes plot(s)\n")
  if (complexity_plot) {
    all_plots$detected_plot <- plot_detected_curves(cpm.mat, FALSE)
  }
  if (complexity_plot_2) {
    all_plots$detected_plot_2 <- plot_detected_curves(cpm.mat, TRUE)
  }

  # ERCC raw reads scatterplot
  if (ERCC_mix) {
    cat("Making ERCC scatter plots\n")
    if (sum(sample.metadata$ERCC == 1) > 0) {
      all_plots$ercc_mix1_plot <- plot_ERCC_observed_scatter(arnaq.run$ERCC.data.cpm,
                                                             sample.metadata, 1)
    } else {
      all_plots$ercc_mix1_plot <- NULL
    }
    if (sum(sample.metadata$ERCC == 2) > 0) {
      all_plots$ercc_mix2_plot <- plot_ERCC_observed_scatter(arnaq.run$ERCC.data.cpm,
                                                             sample.metadata, 2)
    } else {
      all_plots$ercc_mix2_plot <- NULL
    }
  }

  # ERCC fold-change estimates
  if (ERCC_fc) {
    cat("Making ERCC f/c plots\n")
    all_plots$ercc_fc_ps <- plot_all.ERCC.pairs(arnaq.run$ERCC.data.cpm, ERCC.pairs, ERCC.combined)
  }

  # Biotype barchart and pie chart
  if (do.biotype) {
    cat("Making biotypes plots\n")
    all_plots$biotype_ps <- plot_counts_barplot(biotype.counts$reads, 
                                                "Read biotype counts by sample",
                                                max.samples.per.page)
    all_plots$biotype_pie_plot <- plot_pie_biotype(
      biotype.counts$genes,
      "Genes with non-zero counts in at least one sample"
    )
  }

  # Scatter pairs
  if (scatterplots) {
    cat("Making read scatterplots\n")
    all_plots$read_scatterplots_ps <- plot_all_scatter_pairs(norm.data, count.data,
      scatter.pairs, log = TRUE
    )
  }

  assign("all_plots", all_plots, 1)

  # SVG export if required
  if (svg.export) {
    export_all_svg(paste0(plot_directory, "/svg"), qc.name)
  }

  # Save this so you can rerun report generation more easily.
  arnaq.plot.flags <<- list(vsn_normalised=(normalised == "VSN"),
                            linear_normalised=(normalised == "linear"),
                            complexity_plot=complexity_plot, complexity_plot_2=complexity_plot_2,
                            ERCC_show=ERCC_show, ERCC_mix=ERCC_mix, ERCC_fc=ERCC_fc,
                            do.biotype=do.biotype, scatterplots=scatterplots,
                            featureCount.metrics=featureCount.metrics,
                            picard.metrics=picard.metrics)

  arnaq.create.report(rmarkdown_template, qc.name, genome.file, arnaq.plot.flags)
}

#' Creates an ARNAQ report from the premade plots.
#' 
#' This creates an html report using the previously generated plots and an appropriate report 
#' template. This is called automatically by `\link{arnaq}`, but can be run separately if you have
#' altered any of the generated plots to suit your needs.
#' 
#' The `arnaq.plot.flags` and `all_plots` objects must exist.
#' 
#' @param rmarkdown_template The path to the report template.
#' @param qc.name The full name of the report (including any `model.name` appended). Any underscores
#' will be converted to spaces.
#' @param genome.file A character object that describes the genome used. Note this is not used
#' by any data processing, and is used purely so the report includes the genome/version for 
#' completeness. The default uses whatever `\link{arnaq}` found in the `resources.yml` file.
#' @param flags The vector of flags that defines which sections are included in the final report.
#' A suitable vector is automatically created during the main `/link{arnaq}` run, under the name
#' `arnaq.plot.flags`, which is the default.
#' @seealso `\link{arnaq}`
#' 
#' @export
arnaq.create.report <- function(rmarkdown_template = arnaq.report.template, qc.name = qc.name, 
                                genome.file = genome.file, flags = arnaq.plot.flags) {
  if (rmarkdown::pandoc_available()) {
    qc.name.2 <- gsub("_", " ", qc.name)
    rmdfile <- preprocess.rmd(rmarkdown_template, qc.name.2, genome.file, flags)
    render.rmd(rmdfile, qc.name)
    file.remove(rmdfile) # Cleanup after rendering
  } else {
    stop("Pandoc not available!")
  }
}