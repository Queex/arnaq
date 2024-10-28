###################################
# Processing

# Makes gene masks. 
make.gene.masks <- function(count.data, species.gtf, sample.mask = TRUE, ERCC = TRUE) {
  out <- list()

  if (!is.null(species.gtf)) {
    out$Genes_Z <- row.names(count.data) %in% species.gtf$gene_id
  } else {
    tmp1 <- grepl("ERCC-", rownames(count.data), invert=TRUE) # Remove any ERCC
    tmp2 <- grepl("^__", rownames(count.data), invert=TRUE) # Remove any htseq.count summary lines
    out$Genes_Z <- tmp1 & tmp2
  }
  names(out$Genes_Z) <- row.names(count.data)
  cat("Total lines in count table:", length(out$Genes_Z), "\n")
  cat("Total after removing non-gene lines:", sum(out$Genes_Z), "\n")

  # Remove zero lines
  out$Genes <- rowSums(count.data[, sample.mask]) > 0 & out$Genes_Z
  cat(paste("Total genes with >0 reads:", sum(out$Genes), "\n"))

  if (length(out$Genes) == 0) {
    stop("No genes with non-zero counts!")
  }

  if (!is.null(species.gtf)) {
    tmp.gtf <- species.gtf[match(rownames(count.data), species.gtf$gene_id), ]

    out <- make.gene.masks.helper(out, tmp.gtf, "protein_coding")
    out <- make.gene.masks.helper(out, tmp.gtf, "miRNA")
    out <- make.gene.masks.helper(out, tmp.gtf, "lncRNA")
    out <- make.gene.masks.helper(out, tmp.gtf, "sRNA") # tRNA
    out <- make.gene.masks.helper(out, tmp.gtf, "snoRNA")
    out <- make.gene.masks.helper(out, tmp.gtf, "snRNA")
    out <- make.gene.masks.helper(out, tmp.gtf, "scaRNA")
    out$ncRNA_Z <-
      out$miRNA_Z | out$lncRNA_Z | out$sRNA_Z | out$snoRNA_Z | out$snRNA_Z | out$scaRNA_Z
    out$ncRNA <- out$miRNA | out$lncRNA | out$sRNA | out$snoRNA | out$snRNA | out$scaRNA
  }

  if (ERCC) {
    out$ERCC <- grep("ERCC-", rownames(count.data))
    cat("ERCC spike-ins:", length(out$ERCC), "\n")
  }

  return(out)
}

# Helper function for creating gene masks based on biotype
make.gene.masks.helper <- function(masks, tmp.gtf, biotype, sample.mask) {
  masks[[paste0(biotype, "_Z")]] <- tmp.gtf$gene_biotype==biotype
  masks[[paste0(biotype)]] <- masks[[paste0(biotype, "_Z")]] & masks$Genes
  return(masks)
}

# Creates DESeq2 object and performs normalisation, if requested
apply.DESeq2 <- function(masked.data, masked.metadata, form) {

  cat("Creating DESeq2 object\n")
  out <- DESeq2::DESeqDataSetFromMatrix(masked.data, masked.metadata, form)
  out <- DESeq2::estimateSizeFactors(out)
  out <- DESeq2::estimateDispersions(out)
  out <- DESeq2::varianceStabilizingTransformation(out)
  #vsd.data <- as.data.frame(out@assays@data@listData)
  vsd.data <- SummarizedExperiment::assay(out)

  return(list(out, vsd.data))
}