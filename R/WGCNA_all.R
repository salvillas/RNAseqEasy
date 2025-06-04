#' WGCNA Modules Analysis
#'
#' This function performs WGCNA analysis to identify modules of co-expressed genes and their association with traits.
#'
#' @param CurrentDir Directory where results will be saved.
#' @param samplesDir Directory containing sample folders with quant.sf files.
#' @param sample_table Data frame with sample metadata.
#' @param Selection Named vector indicating samples to include/exclude.
#' @param DEGs Vector of differentially expressed genes.
#' @param Variables Character vector of column names in sample_table used for grouping.
#' @param Filter Character vector of samples to exclude.
#' @param Power Soft-thresholding power for WGCNA.
#' @param Name Character. Base name for output files.
#' @param Colors_plot Named vector of colors for groups.
#' @param NumberCol Integer. Number of columns for facet_wrap in ggplot.
#' @param Plot_Ancestors Logical. Whether to plot ancestor GO terms.
#'
#' @return A list with WGCNA results and plots.
#' @export
WGCNA_Modules <- function(CurrentDir, samplesDir, sample_table, Selection = "None", DEGs, Variables,
                          Filter = "None", Power, Name, Colors_plot = NULL, NumberCol = 1, Plot_Ancestors = FALSE) {
  # Subset samples
  sample_table_some <- get_sample_subset(sample_table, Include = names(Selection[Selection == "YES"]),
                                         Exclude = names(Selection[Selection == "NO"]))

  # Load expression data
  txi_some <- load_tximport_data(samplesDir, sample_table_some, tx2gene)
  TPM_all <- as.data.frame(txi_some$abundance)
  colnames(TPM_all) <- sapply(colnames(TPM_all), f1)

  # Filter genes
  TPM_filt <- TPM_all[rownames(TPM_all) %in% DEGs, ]

  # Prepare traits and expression matrices
  datTraits <- prepare_traits(sample_table_some, Variables)
  datExpr <- t(TPM_filt)

  # Quality check
  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK) {
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  }

  # Sample clustering
  plot_sample_clustering(datExpr, CurrentDir)

  # Soft-thresholding
  sft <- pickSoftThreshold(datExpr, powerVector = c(1:10, seq(12, 40, 2)), verbose = 5)
  plot_soft_threshold(sft, CurrentDir)

  if (missing(Power)) stop("Choose a power based on softthresholdingpower.pdf")

  # Network construction
  net <- blockwiseModules(datExpr, power = Power, TOMType = "signed", minModuleSize = 30,
                          reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE,
                          pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "Salva", verbose = 3)

  # Save module info and plots
  save_module_info(net, datExpr, datTraits, CurrentDir)
  plot_module_trait_relationships(net, datExpr, datTraits, CurrentDir)
  plot_module_summaries(net, TPM_filt, sample_table, Variables, Colors_plot, NumberCol, Name, CurrentDir)

  # GO enrichment
  topGO_results <- run_topGO_for_modules(net, geneID2GO_Mpo, CurrentDir, Plot_Ancestors)

  return(topGO_results)
}

#' Prepare traits data for WGCNA
#'
#' This function converts categorical variables into binary traits for WGCNA analysis.
#'
#' @param sample_table Data frame with sample metadata.
#' @param Variables Character vector of column names in sample_table used for grouping.
#'
#' @return A data frame with binary trait columns.
#' @export
prepare_traits <- function(sample_table, Variables) {
  datTraits <- sample_table
  elements <- "Name"
  for (Var in Variables) {
    for (Element in unique(sample_table[[Var]])) {
      datTraits[[Element]] <- ifelse(datTraits[[Var]] == Element, 1, 0)
      elements <- c(elements, Element)
    }
  }
  datTraits <- datTraits[, elements, drop = FALSE]
  rownames(datTraits) <- NULL
  return(datTraits)
}


#' Plot sample clustering dendrogram
#'
#' This function generates a hierarchical clustering dendrogram of samples.
#'
#' @param datExpr Expression matrix (samples x genes).
#' @param CurrentDir Directory where the plot will be saved.
#'
#' @return No return value. A PDF file is saved.
#' @export
plot_sample_clustering <- function(datExpr, CurrentDir) {
  sampleTree <- hclust(dist(datExpr), method = "average")
  pdf(file = file.path(CurrentDir, "sampleClustering.pdf"), width = 12, height = 9)
  par(cex = 0.6)
  par(mar = c(0, 4, 2, 0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  abline(h = 15, col = "red")
  dev.off()
}


#' Plot soft-thresholding power selection
#'
#' This function plots the scale-free topology fit index for various soft-thresholding powers.
#'
#' @param sft Output from WGCNA::pickSoftThreshold().
#' @param CurrentDir Directory where the plot will be saved.
#'
file is saved.
#' @export
plot_soft_threshold <- function(sft, CurrentDir) {
  cex1 <- 0.9
  pdf(file = file.path(CurrentDir, "softthresholdingpower.pdf"), width = 5, height = 5)
  plot(sft$fitIndices[, 1],
       -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit, signed R^2",
       type = "n", main = "Scale independence")
  text(sft$fitIndices[, 1],
       -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers, cex = cex1, col = "red")
  abline(h = 0.90, col = "red")
  dev.off()
}


#' Save module information to CSV
#'
#' This function saves the gene-module assignments and calculates module-trait correlations.
#'
#' @param net WGCNA network object.
#' @param datExpr Expression matrix (samples x genes).
#' @param datTraits Trait matrix (samples x traits).
#' @param CurrentDir Directory where the output file will be saved.
#'
#' @return No return value. A CSV file is written.
#' @export
save_module_info <- function(net, datExpr, datTraits, CurrentDir) {
  moduleColors <- labels2colors(net$colors)
  MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs <- orderMEs(MEs0)
  moduleTraitCor <- cor(MEs, datTraits, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

  geneInfo <- data.frame(
    Genes = colnames(datExpr),
    moduleColor = moduleColors
  )
  geneInfo <- geneInfo[order(geneInfo$moduleColor), ]
  write.csv(geneInfo, file = file.path(CurrentDir, "geneInfo.csv"), row.names = FALSE)
}


#' Plot module-trait relationships
#'
#' This function generates a heatmap of correlations between module eigengenes and traits.
#'
#' @param net WGCNA network object.
#' @param datExpr Expression matrix (samples x genes).
#' @param datTraits Trait matrix (samples x traits).
#' @param CurrentDir Directory where the plot will be saved.
#'
#' @return No return value. A PDF file is saved.
#' @export
plot_module_trait_relationships <- function(net, datExpr, datTraits, CurrentDir) {
  moduleColors <- labels2colors(net$colors)
  MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs <- orderMEs(MEs0)
  moduleTraitCor <- cor(MEs, datTraits, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) <- dim(moduleTraitCor)

  pdf(file = file.path(CurrentDir, "moduleTraitRelationships.pdf"), width = 12, height = 9)
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1, 1),
                 main = "Module-trait relationships")
  dev.off()
}


#' Plot module expression summaries
#'
#' This function generates a faceted plot showing the average expression of each module across conditions.
#'
#' @param net WGCNA network object.
#' @param TPM_filt Filtered TPM expression matrix (genes x samples).
#' @param sample_table Data frame with sample metadata.
#' @param Variables Character vector of grouping variables.
#' @param Colors_plot Named vector of colors for groups.
#' @param NumberCol Number of columns in the facet layout.
#' @param Name Base name for output file.
#' @param CurrentDir Directory where the plot will be saved.
#'
#' @return No return value. A PDF file is saved.
#' @export
plot_module_summaries <- function(net, TPM_filt, sample_table, Variables, Colors_plot, NumberCol, Name, CurrentDir) {
  moduleColors <- labels2colors(net$colors)
  geneInfo <- read.csv(file.path(CurrentDir, "geneInfo.csv"))
  Label_Module <- character(0)
  Module_summary <- matrix(0, nrow = 0, ncol = ncol(TPM_filt))
  colnames(Module_summary) <- colnames(TPM_filt)

  for (color in unique(geneInfo$moduleColor)) {
    genes <- geneInfo$Genes[geneInfo$moduleColor == color]
    counts <- TPM_filt[rownames(TPM_filt) %in% genes, , drop = FALSE]
    norm_counts <- normalize_df(counts)
    summary <- colMeans(norm_counts)
    Module_summary <- rbind(Module_summary, summary)
    rownames(Module_summary)[nrow(Module_summary)] <- capitalize(color)
    write.table(counts, file = file.path(CurrentDir, paste0(color, "_genes.txt")),
                quote = FALSE, sep = "\t", row.names = TRUE)
    Label_Module[capitalize(color)] <- paste(capitalize(color), "(", nrow(counts), " genes)")
  }

  Module_summary_df <- as.data.frame(Module_summary)
  Module_summary_df$Module <- rownames(Module_summary_df)
  Module_summary_long <- tidyr::pivot_longer(Module_summary_df, -Module, names_to = "Sample", values_to = "value")
  Module_summary_long <- merge(Module_summary_long, sample_table, by.x = "Sample", by.y = "Name")

  if (is.null(Colors_plot)) {
    unique_groups <- unique(sample_table[[Variables[1]]])
    palette <- RColorBrewer::brewer.pal(length(unique_groups), "Set2")
    Colors_plot <- setNames(palette, unique_groups)
  }

  p <- ggplot(Module_summary_long, aes_string(x = Variables[2], y = "value", group = Variables[1])) +
    geom_point(aes_string(color = Variables[1]), size = 3, alpha = 0.5) +
    geom_line(stat = "summary", fun = mean, aes_string(color = Variables[1]), size = 1) +
    scale_color_manual(values = Colors_plot) +
    facet_wrap(~Module, scales = "free_y", ncol = NumberCol, labeller = labeller(Module = Label_Module)) +
    labs(title = paste("WGCNA Modules Summary", Name), y = "Normalized expression") +
    theme_bw() +
    theme(strip.text = element_text(face = "bold"), legend.position = "bottom")

  ggsave(file.path(CurrentDir, paste0("Summary_Modules_", Name, ".pdf")), plot = p, width = 10, height = 6)
}


#' Load expression data for a module
#'
#' @param color Module color.
#' @param TPM_filt Filtered TPM expression matrix.
#' @param geneInfo Data frame with gene-module assignments.
#'
#' @return A data frame with expression values for genes in the module.
#' @export
load_module_counts <- function(color, TPM_filt, geneInfo) {
  genes <- geneInfo$Genes[geneInfo$moduleColor == color]
  TPM_filt[rownames(TPM_filt) %in% genes, , drop = FALSE]
}


#' Normalize expression matrix by gene
#'
#' This function scales each gene's expression values between 0 and 1.
#'
#' @param df A numeric matrix or data frame.
#'
#' @return A normalized matrix.
#' @export
normalize_df <- function(df) {
  t(apply(df, 1, function(x) (x - min(x)) / (max(x) - min(x))))
}


#' Run topGO enrichment analysis for each WGCNA module
#'
#' This function performs GO enrichment analysis using topGO for each module identified by WGCNA.
#'
#' @param net WGCNA network object.
#' @param geneID2GO Named list mapping gene IDs to GO terms.
#' @param CurrentDir Directory where results will be saved.
#' @param Plot_Ancestors Logical. Whether to generate ancestor plots.
#'
#' @return A named list with GO enrichment results per module.
#' @export
run_topGO_for_modules <- function(net, geneID2GO, CurrentDir, Plot_Ancestors = FALSE) {
  geneInfo <- read.csv(file.path(CurrentDir, "geneInfo.csv"))
  results_list <- list()

  for (color in unique(geneInfo$moduleColor)) {
    genes <- geneInfo$Genes[geneInfo$moduleColor == color]
    geneUniverse <- names(geneID2GO)
    geneList <- factor(as.integer(geneUniverse %in% genes))
    names(geneList) <- geneUniverse

    GOdata <- new("topGOdata",
                  ontology = "BP",
                  allGenes = geneList,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO)

    resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
    allRes <- GenTable(GOdata, weight01 = resultFisher, topNodes = 100)
    allRes <- allRes[as.numeric(allRes$weight01) <= 0.05, ]

    if (nrow(allRes) > 0) {
      allRes$Enrichment <- (as.numeric(allRes$Significant) / as.numeric(allRes$Annotated)) /
        (sum(geneList == 1) / length(geneList))
      allRes <- allRes[allRes$Enrichment > 1, ]
      write.table(allRes, file = file.path(CurrentDir, paste0("topGO_", color, "_BP.txt")),
                  sep = "\t", quote = FALSE, row.names = FALSE)

      # Optional: plot enrichment
      p <- ggplot(allRes, aes(x = Enrichment, y = reorder(Term, Enrichment), size = Significant, color = as.numeric(weight01))) +
        geom_point() +
        scale_color_gradient(low = "#253494", high = "#edf8b1") +
        labs(title = paste("GO Enrichment -", color), x = "Enrichment", y = NULL) +
        theme_minimal()
      ggsave(file.path(CurrentDir, paste0("topGO_", color, "_plot.pdf")), plot = p, width = 8, height = 6)

      results_list[[color]] <- allRes
    }
  }

  return(results_list)
}

