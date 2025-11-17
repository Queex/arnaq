#' Creates a list of PCA plots (by sample name) up to a certain depth.
#' 
#' The list will contain PCA plots of component 1 against each of the other components up to the 
#' supplied depth. Each sample is represented by their name.
#' 
#' @param PC.norm Normalised PC matrix.
#' @param PC.vars Vector of variability explained by the PCs.
#' @param depth The depth up to which to generate PCA plots.
#' @returns A list of PCA plots up to the requested depth.
#' @seealso `\link{plot_PCA_ids}`
#' @seealso `\link{plot_PCA_group_set}`
#' 
#' @export
plot_PCA_ids_set <- function(PC.norm, PC.vars, depth = 2) {
  out <- list()
  i <- 1
  for (d in 2:depth) {
    out[[i]] <- plot_PCA_ids(PC.norm, PC.vars, c(1, d))
    i <- i + 1
  }
  return(out)
}

#' Creates a PCA plot (by sample name).
#' 
#' Each sample is represented by their name.
#' 
#' @param PC.norm Normalised PC matrix.
#' @param PC.vars Vector of variability explained by the PCs.
#' @param components The two components to include.
#' @returns A PCA plot.
#' @seealso `\link{plot_PCA_ids_set}`
#' @seealso `\link{plot_PCA_groups}`
#' 
#' @export
plot_PCA_ids <- function(PC.norm, PC.vars, components = 1:2) {
  PC.norm$labs <- rownames(PC.norm)
  axis.names <- paste0("PC", components)

  the.plot <- ggplot2::ggplot(PC.norm, 
                              ggplot2::aes(x = .data[[axis.names[1]]],
                                           y = .data[[axis.names[2]]])) +
    ggplot2::geom_text(ggplot2::aes(label = labs), size = 3) +
    ggplot2::labs(
      x = paste0(axis.names[1], " (", PC.vars[components[1]], "%)"),
      y = paste0(axis.names[2], " (", PC.vars[components[2]], "%)"),
      title = "PCA"
    )

  return(the.plot)
}

#' Creates a list of lists of PCA plots (coloured by group membership, for provided groups) up to a 
#' certain depth.
#' 
#' The top-level list will contain one entry per group examined. The sublists will contains PCA 
#' plots of component 1 against each of the other components up to the 
#' supplied depth. Each sample is coloured by group membership of the group for that plot.
#' 
#' @param PC.norm Normalised PC matrix.
#' @param PC.vars Vector of variability explained by the PCs.
#' @param sample.metadata The metadata with group membership for the samples.
#' @param cols A vector naming which of the columns in `sample.metadata` should be used to create
#' PCA plots.
#' @param depth The depth up to which to generate PCA plots. 
#' @returns A list of PCA plots up to the requested depth.
#' @seealso `\link{plot_PCA_group}`
#' @seealso `\link{plot_PCA_group_set}`
#' @seealso `\link{plot_PCA_ids_set}`
#' 
#' @export
plot_all_PCA_groups <- function(PC.norm, PC.vars, sample.metadata, cols, depth = 2) {
  out <- list()
  for (cc in cols) {
    gg <- colnames(sample.metadata)[cc]
    out[[gg]] <- plot_PCA_group_set(PC.norm, PC.vars, sample.metadata[, cc], gg, depth = depth)
  }
  return(out)
}

#' Creates a list of PCA plots (coloured by group membership) up to a certain depth.
#' 
#' The list will contain PCA plots of component 1 against each of the other components up to the 
#' supplied depth. Each sample is coloured by group membership of a specific group.
#' 
#' @param PC.norm Normalised PC matrix.
#' @param PC.vars Vector of variability explained by the PCs.
#' @param groups Vector of group membership for the samples.
#' @param groupName The name for this group, displayed in the legend.
#' @param depth The depth up to which to generate PCA plots. 
#' @returns A list of PCA plots up to the requested depth.
#' @seealso `\link{plot_PCA_group}`
#' @seealso `\link{plot_all_PCA_groups}`
#' @seealso `\link{plot_PCA_ids_set}`
#' 
#' @export
plot_PCA_group_set <- function(PC.norm, PC.vars, groups, groupName, depth = 2) {
  out <- list()
  i <- 1
  for (d in 2:depth) {
    out[[i]] <- plot_PCA_group(PC.norm, PC.vars, groups, groupName, c(1, d))
    i <- i + 1
  }
  return(out)
}

#' Creates a PCA plot (coloured by group).
#' 
#' Each sample is coloured by group membership.
#' 
#' @param PC.norm Normalised PC matrix.
#' @param PC.vars Vector of variability explained by the PCs.
#' @param groups Vector of group membership for the samples.
#' @param groupName The name for this group, displayed in the legend.
#' @param components The two components to include.
#' @returns A PCA plot.
#' @seealso `\link{plot_PCA_ids}`
#' @seealso `\link{plot_PCA_group_set}`
#' @seealso `\link{plot_all.PCA.group.set}`
#' 
#' @export
plot_PCA_group <- function(PC.norm, PC.vars, groups, groupName, components = 1:2) {
  PC.norm[[groupName]] <- groups
  axis.names <- paste0("PC", components)

  the.plot <- ggplot2::ggplot(PC.norm, ggplot2::aes(
    x = .data[[axis.names[1]]], y = .data[[axis.names[2]]],
    colour = .data[[groupName]], shape = .data[[groupName]]
  )) +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = paste0(axis.names[1], " (", PC.vars[components[1]], "%)"),
      y = paste0(axis.names[2], " (", PC.vars[components[2]], "%)"),
      title = paste("PCA for group", groupName)
    )

  tmp2 <- length(unique(groups))
  if (tmp2 > 6) {
    tmp <- rep(c(16, 17, 15, 3, 7, 8), times = ceiling(tmp2 / 6))[1:tmp2]
    the.plot <- the.plot +
      ggplot2::scale_shape_manual(values = tmp)
  }

  return(the.plot)
}

# Helper function to calculate PCA, returning a list of useful values.
calc.pca <- function(data) {
  cat("Calculating PCA\n")
  pca.norm <- stats::prcomp(t(data))
  PC.norm <- data.frame(pca.norm$x)

  var.cap <- summary(pca.norm)
  pc_var.cap <- round(var.cap$importance[2, ] * 100, 1)
  cum_var.cap <- round(var.cap$importance[3, ] * 100, 1)

  return(list(PC = PC.norm, Vars = pc_var.cap, CumulVars = cum_var.cap))
}