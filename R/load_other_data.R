# Reads ERCC concentration file
read.ERCC.table <- function(source_vec) {
  cat("Reading ERCC table\n")
  ERCC.table <<- fallback_load_table(source_vec, sep = "\t", header = TRUE)
}