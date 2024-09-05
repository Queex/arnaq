# Helper function for saving out count tables sensibly.
save.counts <- function(data, filename) {
  cat(paste("Writing count data to table:", filename, "\n"))
  write.table(data, filename, sep = "\t", quote = FALSE, row.names = TRUE)
}

# Convenience function to make sure a directory exists to write to.
ensure.directory <- function(out.directory) {
  # Ensure target directory exists
  dir.create(out.directory, showWarnings = FALSE)
}

# Ensure no sinks are left behind when starting run
clear.sinks <- function() {
  while (sink.number()>0) {
    sink()
  }
}