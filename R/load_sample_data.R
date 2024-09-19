# Loads count data
read.count.data <- function(sources_vec) {
  cat("Reading count data\n")
  df <- fallback_load_table(sources_vec, row.names = 1, header = TRUE)
  cat("Sample Names: ")
  cat(colnames(df))
  cat("\n")

  cat("Data dimensions: ")
  cat(dim(df))
  cat("\n")

  return(df)
}

# Loads the assignment summary file
read.read.summary <- function(sources_vec, count.data, gene.masks, sample.metadata) {
  cat("Reading assignment summary\n")
  read.summary <- fallback_load_table(sources_vec, header = TRUE, row.names = 1)
  if (!("Genes" %in% rownames(read.summary))) { # Only do this fix if not already done
    row.names(read.summary) <- gsub("Unassigned_", "", row.names(read.summary))
    genes <- apply(count.data[gene.masks$Genes, ], 2, sum)
    ERCC <- apply(count.data[gene.masks$ERCC, ], 2, sum)
    non_gene <- read.summary["Assigned", ] - genes - ERCC
    read.summary <- rbind(Genes = genes, Non_Gene = non_gene, ERCC = ERCC, read.summary)
    read.summary <- read.summary[-4, ] # Should only ever be Assigned row
    read.summary <- read.summary[rowSums(read.summary) > 0, ] # Remove 0 rows
  }

  if (!identical(colnames(read.summary), colnames(count.data))) {
    if (!identical(colnames(read.summary), sample.metadata$Name)) {
      cat(colnames(read.summary))
      cat("\n")
      cat(colnames(count.data))
      cat("\n")
      stop("Error in samples for read summary")
    }
  }
  read.summary
}

# Loads the collated Picard metrics files
read.picard.metrics <- function(picard_sources_vec, dups_sources_vec, sample.names = NULL) {
  cat("Reading Picard metrics\n")
  tmp <- fallback_load_table(picard_sources_vec, header = TRUE, sep = "\t", row.names = 1)
  if (is.data.frame(tmp)) { # Only do fix up if reading from file
    base.stats <- tmp[, c(
      "PF_BASES", "PF_ALIGNED_BASES", "CODING_BASES",
      "UTR_BASES", "INTRONIC_BASES", "INTERGENIC_BASES"
    )]
    pct.stats <- tmp[, c(
      "MEDIAN_CV_COVERAGE", "MEDIAN_5PRIME_BIAS",
      "MEDIAN_3PRIME_BIAS", "MEDIAN_5PRIME_TO_3PRIME_BIAS"
    )]
    cat("Reading duplication rates\n")
    tmp <- fallback_load_table(dups_sources_vec, header = TRUE, sep = "\t", row.names = 1,
                               nrow=nrow(base.stats)) # This remotes the histogram portion
    pct.stats <- cbind(pct.stats, tmp$PERCENT_DUPLICATION)
    colnames(base.stats) <- c(
      "All", "Aligned", "Coding", "UTR", "Intronic",
      "Intergenic"
    )
    colnames(pct.stats) <- c(
      "Median CV Coverage", "Median 5' Bias",
      "Median 3' Bias", "Median Bias Ratio", "Duplication Rate"
    )
    if (!is.null(sample.names)) {
      rownames(pct.stats) <- rownames(base.stats) <- sample.names
    }
  }

  return(list(ByBases = base.stats, ByOther = pct.stats))
}