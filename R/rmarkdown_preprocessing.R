# Preprocess the .Rmd template.
preprocess.rmd <- function(in.filename, project.name, genome.file, plot.flags) {
  cat("Preprocessing rmd file\n")
  rmd <- paste(readLines(in.filename), collapse = "\n")

  tags <- list(vsn_normalised = "Normalisation_VSN_section",
    linear_normalised = "Normalisation_linear_section",
    complexity_plot = c("Complexity_plot1_section", "Detected_plot1_section"),
    complexity_plot_2 = c("Complexity_plot2_section", "Detected_plot2_section"),
    ERCC_show = "ERCC_section", ERCC_mix = "ERCC_mix_reads", ERCC_fc = "ERCC_fc_accuracy",
    do.biotype = "Biotype_section", scatterplots = "Scatterplot_section",
    featureCount.metrics = "Barplot_section", picard.metrics = "Violin_section"
  )

  # Remove unneeded sections
  for (i in seq_along(plot.flags)) {
    if (!plot.flags[[i]]) {
      for (to.remove in tags[[names(plot.flags)[i]]]) {
        cat(paste("Removing", to.remove, "\n"))
        rmd <- preprocess.remove(rmd, to.remove)
      }
    }
  }

  # Replace variables
  rmd <- preprocess.sub(rmd, "Project_number", project.name)
  rmd <- preprocess.sub(rmd, "Date", Sys.Date())
  rmd <- preprocess.sub(rmd, "num_samples", ncol(norm.data))
  rmd <- preprocess.sub(rmd, "genome_file", genome.file)
  rmd <- preprocess.sub(rmd, "genes", sum(gene.masks$Genes_Z))
  rmd <- preprocess.sub(rmd, "detected_genes", sum(gene.masks$Genes))

  # Looped sections
  if (plot.flags$featureCount.metrics) {
    rmd <- preprocess.loop(rmd, "Count_barplot_loop", length(all_plots$barplots_ps))
    rmd <- preprocess.loop(rmd, "Prop_barplot_loop", length(all_plots$barplots2_ps))
  }

  rmd <- preprocess.loop( # PCA by groups
    rmd, "PCA_group_loop", length(all_plots$pca_group_pss),
    names(all_plots$pca_group_pss)
  )

  if (plot.flags$ERCC_fc) {
    rmd <- preprocess.loop(
      rmd, "ERCC_fc_loop", length(all_plots$ercc_fc_ps),
      names(all_plots$ercc_fc_ps)
    )
  }

  if (plot.flags$do.biotype) {
    rmd <- preprocess.loop(rmd, "Biotype_loop", length(all_plots$biotype_ps))
  }

  if (plot.flags$scatterplots) {
    rmd <- preprocess.loop(
      rmd, "Scatterplot_loop", length(all_plots$read_scatterplots_ps),
      names(all_plots$read_scatterplots_ps)
    )
  }

  # Clean remainder
  rmd <- preprocess.clean(rmd)

  # Save preprocessed rmd
  out.file <- paste0("QC/", project.name, "_arnaq_tmp.rmd")
  cat(rmd, file = out.file, sep = "")

  return(out.file)
}

# Render the preprocessed .Rmd template
render.rmd <- function(rmdfile = arnaq.report.template, qc.name = qc.name) {
  html.file <- paste0(qc.name, "_arnaq.html")
  pdf.file <- paste0(qc.name, "_arnaq.pdf")
  rmarkdown::render(rmdfile,
    output_format = c("html_document", "pdf_document"),
    output_file = c(html.file, pdf.file), output_dir = "QC"
  )
}

# Checks that the template is the correct version.
check.QC.template <- function(in.filename, expected.template.version) {
  cat("Checking template exists...")
  rmd <- paste(readLines(in.filename), collapse = "\n")

  if (!file.exists(in.filename)) {
    cat("\n")
    stop(paste("Template ", in.filename, "does not exist."))
  }
  if (!grepl("<!--- Generated ARNAQ document --->", rmd, TRUE)) {
    cat("\n")
    stop(paste(in.filename, "does not seem to be a template for ARNAQ."))
  }
  tmp <- paste0("Version: ", expected.template.version)
  if (!grepl(tmp, rmd, TRUE)) {
    cat("\n")
    stop(paste("This template appears to be for a different version of ARNAQ."))
  }
  cat(" done.\n")
}

# Removes a tagged section from the template.
preprocess.remove <- function(rmd, tag) {
  # cat(paste("Removing:",tag,"\n"))
  return(sub(paste0("@@", tag, ".*@@/", tag), "", rmd))
}

# Replaces a tag with a value.
preprocess.sub <- function(rmd, tag, val) {
  # cat(paste("Replacing",tag,"with",val,"\n"))
  return(gsub(paste0("\\$", tag), val, rmd))
}

# Replaces a tagged section with a list of numbered subsections
preprocess.loop <- function(rmd, tag, iter, names = "") {
  # cat(paste("Loop:",tag,"[",paste(names,collapse=","),"]\n"))
  loopy <- sub(paste0(".*@@", tag), "", rmd)
  loopy <- sub(paste0("@@/", tag, ".*"), "", loopy)
  # Remove header if singleton plot
  if (iter == 1) {
    loopy <- sub("#+[^`]*", "", loopy)
  }
  out <- ""
  for (i in 1:iter) {
    tmp <- preprocess.sub(loopy, "i", i)
    tmp <- preprocess.sub(tmp, "name", names[i])
    out <- paste(out, tmp, sep = "\n\n")
  }

  out <- sub(paste0("@@", tag, ".*@@/", tag), out, rmd)
  return(out)
}

# Finds any remaining tags and removes them.
preprocess.clean <- function(rmd) {
  return(gsub("@@[^\n]*\n", "", rmd))
}