# Reads the sample file.
read.samples.metadata <- function(source_vec = "samples.txt", count.data) {
  metadata <- fallback_load_table(source_vec, header = TRUE)
  if (ncol(metadata) < 4) {
    stop("`samples.txt` has too few columns.")
  }
  for (i in 4:ncol(metadata)){
    metadata[, i] <- as.factor(metadata[, i])
  }

  if (!identical(colnames(count.data), as.character(metadata$Display))) {
    # Only sort if not using cached data
    metadata <- metadata[match(as.character(colnames(count.data)), as.character(metadata$Name)), ]
    rownames(metadata) <- metadata$Display
  } else {
    # Now test that the sorting has worked and no rows are lost
    if (!identical(as.character(metadata$Display), as.character(colnames(count.data)))) {
      cat(as.character(metadata$Name))
      cat("\n")
      cat(as.character(colnames(count.data)))
      cat("\n")
      stop("Error reconciling sample names!")
    }
  }

  cat("Updated names: ")
  cat(as.character(metadata$Display))
  cat("\n")
  return(metadata)
}

# Reads the resources file.
read.resources.file <- function(filename = "resources.yml") {
  cat("Reading resources file\n")
  out <- list()
  res.file <- file(filename, open = "r")
  for (line in readLines(res.file)) {
    tmp <- unlist(strsplit(line, ":"))
    out[[trimws(tmp[1])]] <- trimws(tmp[2])
  }
  cat("Settings found -\n")
  for (n in names(out)) {
    cat(paste0(n, ": ", out[[n]]), "\n")
  }
  close(res.file)
  return(out)
}