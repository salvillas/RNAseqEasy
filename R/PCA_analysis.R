#' Calculate PCA from a DESeq2 object
#'
#' This function performs variance-stabilizing transformation (VST),
#' selects the top variable genes, and computes PCA.
#'
#' @param dds A DESeqDataSet object.
#' @param ntop Integer. Number of top variable genes to use. Default is 500.
#' @return A list with PCA results, percent variance, and the transformed object.
calculate_pca <- function(dds, ntop = 500) {
  vsd <- DESeq2::vst(dds, blind = FALSE)
  rv <- matrixStats::rowVars(assay(vsd))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(vsd)[select, ]))
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  list(pca = pca, percentVar = percentVar, vsd = vsd)
}

#' Plot PCA from PCA results
#'
#' This function generates a PCA plot using ggplot2.
#'
#' @param pca_result Output from `calculate_pca()`.
#' @param variables Character vector of column names in colData(vsd) to use for grouping.
#' @param colors Named vector of colors for groups.
#' @param components Integer vector of length 2 indicating which PCs (1 to 4) to plot
#' @param name Character. Name used for saving the plot.
#' @param output_dir Directory where the plot will be saved.
#' @param width Width units for the plot. Default to 10.
#' @param height Height units for the plot. Default to 8.
#' @param units Units to set plot dimensions. Default to "in" (inches)
#' @return A ggplot object.
plot_pca <- function(pca_result, variables, colors, components = c(1, 2), name = "PCA",
                     output_dir = ".", width = 10, height = 8, units = "in") {
  vsd <- pca_result$vsd
  pca <- pca_result$pca
  percentVar <- pca_result$percentVar
  if (!all(variables %in% names(colData(vsd)))) {
    stop("The argument 'variables' must match columns in colData(vsd).")
  }

  intgroup.df <- as.data.frame(colData(vsd)[, variables, drop = FALSE])
  group <- if (length(variables) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  } else {
    colData(vsd)[[variables]]
  }

  d <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3],
    PC4 = pca$x[, 4],
    group = group,
    intgroup.df,
    name = colnames(vsd)
  )

  attr(d, "percentVar") <- percentVar[1:2]

  p <- ggplot(d, aes_string(x = paste0("PC", components[1]), y = paste0("PC", components[2]))) +
    geom_point(size = 3, aes_string(color = variables[2], shape = variables[1])) +
    scale_color_manual(values = colors) +
    xlab(paste0("PC", components[1], ": ", round(percentVar[components[1]] * 100), "% variance")) +
    ylab(paste0("PC", components[2], ": ", round(percentVar[components[2]] * 100), "% variance")) +
    coord_fixed() +
    theme_light()

  ggsave(file.path(output_dir, paste0("PCA_plot_", name, ".png")), p, width = width, height = height, units = units, limitsize = FALSE)

  return(p)
}

#' Run full PCA analysis and save plot
#'
#' This function runs the full PCA pipeline: transformation, PCA computation, and plotting.
#'
#' @param dds A DESeqDataSet object.
#' @param variables Character vector of column names in colData(dds) to use for grouping.
#' @param colors Named vector of colors for groups.
#' @param components Integer vector of length 2 indicating which PCs (1 to 4) to plot
#' @param name Character. Name used for saving the plot.
#' @param output_dir Directory where the plot will be saved.
#' @param ntop Integer. Number of top variable genes to use. Default is 500.
#' @param width Width units for the plot. Default to 10.
#' @param height Height units for the plot. Default to 8.
#' @param units Units to set plot dimensions. Default to "in" (inches)
#' @return A ggplot object.
run_pca_analysis <- function(dds, variables, colors, components = c(1, 2), name = "PCA", output_dir = ".", ntop = 500,
                             width = 10, height = 8, units = "in") {
  pca_result <- calculate_pca(dds, ntop = ntop)
  plot <- plot_pca(pca_result, variables, colors, components, name, output_dir, width = width, height = height, units = units)
  return(plot)
}
