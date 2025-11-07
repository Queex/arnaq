# Reads the biotype mapping table
read.biotype.conversion.table <- function(sources_vec) {
  cat("Reading biotype conversion table\n")
  tmp <- fallback_load_table(sources_vec, header = FALSE, sep = "\t", quote = "")
  row.names(tmp) <- tmp[, 1]
  tmp
}

# Reads the gtf reference file, including biotypes
read.biotypes <- function(annotation.file) {
  cat("Reading gtf file\n")
  if (exists("species.gtf")) {
    species.gtf
  } else {
    species.gtf <- rtracklayer::readGFF(annotation.file)
    species.gtf <- species.gtf[species.gtf$type == "gene", ]
    species.gtf
  }
}