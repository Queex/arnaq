# Reads the biotype mapping table
read.biotype.conversion.table <- function(sources_vec) {
  cat("Reading biotype conversion table\n")
  tmp <- fallback_load_table(sources_vec, header = FALSE, sep = "\t", quote = "")
  row.names(tmp) <- tmp[, 1]
  tmp
}

# Reads the gtf reference file, including biotypes
read.annotation <- function(annotation.dir, count.data) {
  cat("Reading gtf file\n")
  annotation.file <- paste(annotation.dir, list.files(annotation.dir, pattern = ".*\\.gtf$")[1],
                           sep = "/")
  species.gtf <- rtracklayer::readGFF(annotation.file)
  species.gtf <- species.gtf[species.gtf$type == "gene", ]
  species.gtf <- species.gtf[match(rowNames(count.data), species.gtf$gene_name), ]
  row.names(species.gtf) <- species.gtf$gene_name
  species.gtf
}