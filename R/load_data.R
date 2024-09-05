# Function that tries to load data from a number of different sources, in order.

fallback_load_table <- function(sources_vec, ...) {
  for (source in sources_vec) {
    if (exists(source)) {
      return(get(source)) # if named object exists, use this
    } else if (file.exists(source)) {
      return(read.table(source, ...)) # if file exists, read table
    }
  }
  stop(paste("No objects or files found for fallback load:", paste(sources_vec, collapse=", ")))
}