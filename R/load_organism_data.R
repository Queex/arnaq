# Reads the biotype mapping table
read.biotype.conversion.table <- function(sources_vec) {
  cat("Reading biotype conversion table\n")
  tmp <- fallback_load_table(sources_vec, header = FALSE, sep = "\t", quote = "")
  row.names(tmp) <- tmp[, 1]
  return(tmp)
}

# Reads the gtf reference file, including biotypes
read.biotypes <- function(species, resource.dir) {
  cat("Reading gtf file\n")
  if (exists("species.gtf")) {
    return(species.gtf)
  } else {
    annotation.dir <- paste(resource.dir, species, sep = "/")
    annotation.file <- paste(annotation.dir, list.files(annotation.dir, pattern = ".*\\.gtf$")[1],
                             sep = "/")
    species.gtf <- rtracklayer::readGFF(annotation.file)
    species.gtf <- species.gtf[species.gtf$type == "gene", ]
    return(species.gtf)
  }
}