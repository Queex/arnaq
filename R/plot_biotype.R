#' Create a pie chart showing the distribution of biotypes across observed transcripts.
#' 
#' This is not weighted by the abundance of reads, and only serves as a breakdown of the biotypes 
#' captured in the data set.
#' @param gene_counts Vector of gene counts, named with the biotype.
#' @param main The title to give the plot.
#' @returns A ggplot2 plot of a pie chart.
#'
#' @export
plot_pie_biotype <- function(gene_counts, main) {
  group.names <- names(gene_counts)
  group.names <- paste0(group.names, " (", gene_counts, ")")
  cat.name <- "Biotype"

  tmp <- data.frame(names = group.names, value = as.numeric(gene_counts))

  return(ggplot2::ggplot(
    tmp,
    ggplot2::aes(x = factor(1), y = value, fill = factor(group.names, levels = group.names))
  ) +
    ggplot2::geom_col(width = 1, position = "fill") +
    ggplot2::coord_polar(theta = "y", start = 0) +
    ggplot2::labs(x = "", y = "", fill = cat.name, title = main) +
    ggplot2::theme(
      axis.ticks = ggplot2::element_line(colour = "#00000000"),
      axis.text.y = ggplot2::element_text(colour = "#00000000")
    ))
}

# Turns gene count data into biotype count data. Returns a list with 3 members: read counts,
# gene counts, and reads without collating into sensible categories.
make.biotype.data <- function(count.data, species.gtf, biotype.conversion.table) {
  dat <- as.matrix(count.data)

  species.biotype <- species.gtf[, c("gene_id", "gene_biotype")]
  species.biotype <- species.biotype[match(row.names(count.data), species.biotype$gene_id), ]
  row.names(species.biotype) <- species.biotype[, 1]
  biotype <- unique(biotype.conversion.table[, 2])
  biotype.unconverted <- unique(biotype.conversion.table[, 1])

  gene.type <- species.biotype[row.names(dat), 2] # change row names to biotypes
  gene.type[is.na(gene.type)] <- "unknown"
  out3 <- matrix(NA, ncol = ncol(dat), nrow = length(biotype.unconverted)) # unconverted read counts
  rownames(out3) <- biotype.unconverted
  for (ii in seq_along(biotype.unconverted)) {
    bt <- biotype.unconverted[ii]
    out3[ii, ] <- colSums(dat[gene.type == bt, , drop = FALSE])
  }

  row.names(biotype.conversion.table) <- biotype.conversion.table[, 1]
  gene.type2 <- biotype.conversion.table[gene.type, 2] # convert to reduced types

  out <- matrix(NA, ncol = ncol(dat), nrow = length(biotype)) # for read counts
  rownames(out) <- biotype
  for (ii in seq_along(biotype)) {
    bt <- biotype[ii]
    out[ii, ] <- colSums(dat[gene.type2 == bt, , drop = FALSE])
  }

  colnames(out) <- colnames(out3) <- colnames(dat)

  out <- out[rowSums(out) > 0, ]
  out3 <- out3[rowSums(out3) > 0, ]
  out2 <- table(gene.type2) # for non-zero gene breakdown

  list(reads = out, genes = out2, readsUnconverted = out3)
}