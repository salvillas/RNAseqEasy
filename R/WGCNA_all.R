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
  elements <- "Sample"
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
#' @param output_path Directory where the plot will be saved.
#'
#' @return No return value. A PDF file is saved.
#' @export
plot_sample_clustering <- function(datExpr, output_path) {
  sampleTree <- hclust(dist(datExpr), method = "average")
  pdf(file = file.path(output_path, "sampleClustering.pdf"), width = 12, height = 9)
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
#' @param output_path Directory where the plot will be saved.
#'
#' @export
plot_soft_threshold <- function(sft, output_path) {
  cex1 <- 0.9
  pdf(file = file.path(output_path, "softthresholdingpower.pdf"), width = 5, height = 5)
  plot(sft$fitIndices[, 1],
       -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit, signed R^2",
       type = "n", main = "Scale independence")
  text(sft$fitIndices[, 1],
       -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = c(1:10, seq(12, 40, 2)), cex = cex1, col = "red")
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
#' @param output_path Directory where the output file will be saved.
#' @param Name Character. Base name for output files.

#'
#' @return No return value. A CSV file is written.
#' @export
save_module_info <- function(net, datExpr, datTraits, output_path, Name) {
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
  write.csv(geneInfo, file = file.path(output_path, paste(Name, "geneInfo.csv", sep = "_")), row.names = FALSE)
}


#' Plot module-trait relationships
#'
#' This function generates a heatmap of correlations between module eigengenes and traits.
#'
#' @param net WGCNA network object.
#' @param datExpr Expression matrix (samples x genes).
#' @param datTraits Trait matrix (samples x traits).
#' @param output_path Directory where the plot will be saved.
#'
#' @return No return value. A PDF file is saved.
#' @export
plot_module_trait_relationships <- function(net, datExpr, datTraits, output_path) {
  moduleColors <- labels2colors(net$colors)
  MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs <- orderMEs(MEs0)
  moduleTraitCor <- cor(MEs, datTraits, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) <- dim(moduleTraitCor)

  pdf(file = file.path(output_path, "moduleTraitRelationships.pdf"), width = 12, height = 9)
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
#' @param output_path Directory where the plot will be saved.
#'
#' @return No return value. A PDF file is saved.
#' @export
plot_module_summaries <- function(net, TPM_filt, sample_table, Variables, Colors_plot, NumberCol, Name, output_path) {
  moduleColors <- labels2colors(net$colors)
  geneInfo <- read.csv(file.path(output_path, paste(Name, "geneInfo.csv", sep = "_")))
  Label_Module <- character(0)
  Module_summary <- matrix(0, nrow = 0, ncol = ncol(TPM_filt))
  colnames(Module_summary) <- colnames(TPM_filt)

  for (color in unique(geneInfo$moduleColor)) {
    genes <- geneInfo$Genes[geneInfo$moduleColor == color]
    counts <- TPM_filt[rownames(TPM_filt) %in% genes, , drop = FALSE]
    norm_counts <- normalize_df(counts)
    summary <- colMeans(norm_counts)
    Module_summary <- rbind(Module_summary, summary)
    rownames(Module_summary)[nrow(Module_summary)] <- Hmisc::capitalize(color)
    write.table(counts, file = file.path(output_path, paste0(color, "_genes.txt")),
                quote = FALSE, sep = "\t", row.names = TRUE)
    Label_Module[Hmisc::capitalize(color)] <- paste(Hmisc::capitalize(color), "(", nrow(counts), " genes)")
  }

  Module_summary_df <- as.data.frame(Module_summary)
  Module_summary_df$Module <- rownames(Module_summary_df)
  Module_summary_long <- tidyr::pivot_longer(Module_summary_df, -Module, names_to = "Sample", values_to = "value")
  Module_summary_long <- merge(Module_summary_long, sample_table, by = "Sample")

  if (is.null(Colors_plot)) {
    unique_groups <- unique(sample_table[[Variables[1]]])
    palette <- RColorBrewer::brewer.pal(length(unique_groups), "Set2")
    Colors_plot <- setNames(palette, unique_groups)
  }

  facets_formula <- if (length(Variables) == 3) {
    as.formula(paste0("Module~", Variable[3]))
  } else {
    as.formula("~Module")
  }

  p <- ggplot2::ggplot(Module_summary_long, ggplot2::aes_string(x = Variables[2], y = "value", group = Variables[1])) +
    ggplot2::geom_point(ggplot2::aes_string(color = Variables[1]), size = 3, shape = 20, alpha = 0.5) +
    ggplot2::geom_line(stat = "summary", fun = mean, ggplot2::aes_string(color = Variables[1]), size = 1) +
    ggplot2::scale_color_manual(values = Colors_plot) +
    ggplot2::facet_wrap(facets_formula,
                        scales = "free_y", ncol = NumberCol, labeller = labeller(Module = Label_Module)) +
    ggplot2::labs(title = paste("WGCNA Modules Summary", Name), y = "Normalized expression") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
    ggplot2::theme_bw() +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour= "black"),
                   axis.ticks.x= ggplot2::element_blank(),
                   legend.position = "bottom",
                   legend.title = ggplot2::element_text(face = "bold"),
                   axis.title = ggplot2::element_text(face = "bold"),
                   title = ggplot2::element_text(face= "bold"),
                   strip.text = ggplot2::element_text(face = "bold", hjust = 0))

  ggsave(file.path(output_path, paste0("Summary_Modules_", Name, ".pdf")), plot = p, units = "in",
         width = 8, height = 3 + length(unique(Module_summary_long$Module)))
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
#' @param geneID2GO Named list mapping gene IDs to GO terms.
#' @param Name Character. Base name for '"Name"_geneInfo.csv' previously created with WGCNA.
#' @param input_path Directory where '"Name"_geneInfo.csv' file was saved.
#' @param output_path_topGO Directory where results will be saved.
#' @param ontology GO ontology to use ("BP", "MF", "CC"). Defaults to "BP".
#' @param algorithm Algorithm for topGO test (default: "weight01").
#' @param statistic Statistical test to use (default: "fisher").
#' @param plot_similarity Logical, whether to analyze and visualize GO term similarity.
#' @param Number_GOs Number of top GO term names to plot in the scatterplot. Defaults to 20.
#' @param orgdb OrgDb package name for similarity analysis. Defaults to "org.At.tair.db".
#' @param semdata Optional precomputed semantic data.
#' @param save_GeneNames Logical, whether to save genes represented by GO terms of interest.
#' @param Annotation Two column data frame. First column must include Gene IDs.
#' Second column must include functional annotation.
#' @param Ontologies Character vector of GO terms of interest to search for.
#'
#' @return A named list with GO enrichment results per module.
#' @export
run_topGO_for_modules <- function(geneID2GO, Name,
                                  input_path, output_path_topGO, ontology = "BP",
                                  algorithm = "weight01",
                                  statistic = "fisher",
                                  plot_similarity = TRUE, Number_GOs = 20,
                                  orgdb = "org.At.tair.db", semdata = NULL,
                                  save_GeneNames = FALSE, Annotation, Ontologies = NULL) {
  geneInfo <- read.csv(file.path(input_path, paste(Name, "geneInfo.csv", sep = "_")))
  results_list <- list()

  for (color in unique(geneInfo$moduleColor)) {
    genes <- geneInfo$Genes[geneInfo$moduleColor == color]
    topGO_module <- topGO_All(DEG = genes, geneID2GO = geneID2GO, name = color,
                              output_dir = output_path_topGO, ontology = ontology,
                              algorithm = algorithm, statistic = statistic,
                              plot_similarity = plot_similarity, Number_GOs = Number_GOs,
                              orgdb = orgdb, semdata = semdata, save_GeneNames = save_GeneNames, Ontologies = Ontologies, Annotation = Annotation)
    results_list[[color]] <- topGO_module
  }
  return(results_list)
  }


#' WGCNA Modules Analysis
#'
#' This function performs WGCNA analysis to identify modules of co-expressed genes and their association with traits.
#'
#' @param output_path Directory where results will be saved.
#' @param sampleDir Directory containing sample folders with quant.sf files.
#' @param sample_table Data frame with sample metadata.
#' @param Selection Named vector indicating samples to include/exclude.
#' @param DEGs Vector of differentially expressed genes.
#' @param Variables Character vector of column names in sample_table used for grouping.
#' @param tx2gene A data frame with columns TXNAME and GENEID (transcript-to-gene mapping).
#' @param Filter Character vector of samples to exclude.
#' @param Power Soft-thresholding power for WGCNA.
#' @param Name Character. Base name for output files.
#' @param Colors_plot Named vector of colors for groups.
#' @param NumberCol Integer. Number of columns for facet_wrap in ggplot.
#' @param geneID2GO Named list mapping gene IDs to GO terms.
#' @param ontology GO ontology to use ("BP", "MF", "CC"). Defaults to "BP".
#' @param algorithm Algorithm for topGO test (default: "weight01").
#' @param statistic Statistical test to use (default: "fisher").
#' @param plot_similarity Logical, whether to analyze and visualize GO term similarity.
#' @param Number_GOs Number of top GO term names to plot in the scatterplot. Defaults to 20.
#' @param orgdb OrgDb package name for similarity analysis. Defaults to "org.At.tair.db".
#' @param semdata Optional precomputed semantic data.
#' @param save_GeneNames Logical, whether to save genes represented by GO terms of interest.
#' @param Annotation Two column data frame. First column must include Gene IDs.
#' Second column must include functional annotation.
#' @param Ontologies Character vector of GO terms of interest to search for.
#'
#' @return A list with WGCNA results and plots.
#' @export
WGCNA_Modules <- function(output_path, sampleDir, sample_table, Include = NULL, Exclude = NULL, DEGs, Variables,
                          tx2gene, Filter = NULL, Power, Name, Colors_plot = NULL, NumberCol = 1,
                          geneID2GO, ontology = "BP", algorithm = "weight01", statistic = "fisher",
                          plot_similarity = TRUE, Number_GOs = 20,
                          orgdb = "org.At.tair.db", semdata = NULL,
                          save_GeneNames = save_GeneNames, Annotation, Ontologies = NULL) {
  # Subset samples
  sample_table_some <- get_sample_subset(sample_table, Include = Include,
                                         Exclude = Exclude)
  sample_table_some <- add_sample_path(sampleDir = sampleDir, sample_table = sample_table_some)
  # Load expression data
  txi_some <- load_tximport_data(sampleDir, sample_table_some, tx2gene)
  TPM_all <- as.data.frame(txi_some$abundance)
  # In each sample name there are different things separated by "_" (name, run, etc...). With this function we will get the first one, that is, the real sample name
  f1 <- function(x) stringr::str_split(x, pattern = "_")[[1]][1]
  colnames(TPM_all) <- sapply(colnames(TPM_all), f1)

  # Filter genes
  TPM_filt <- TPM_all[rownames(TPM_all) %in% DEGs, ]

  # Prepare traits and expression matrices
  datTraits <- prepare_traits(sample_table_some, Variables)


  # Filter samples if you want to exclude them
  if (!is.null(Filter)) {
    if (all(Filter %in% colnames(TPM_filt))) {
      TPM_filt <- TPM_filt[,!names(TPM_filt) %in% Filter]
      datTraits <- datTraits[!datTraits$Sample %in% Filter,]
    }
  }

  datExpr <- t(TPM_filt)
  # Quality check
  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK) {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  }

  # Sample clustering
  plot_sample_clustering(datExpr, output_path)

  # Soft-thresholding
  sft <- pickSoftThreshold(datExpr, powerVector = c(1:10, seq(12, 40, 2)), verbose = 5)
  plot_soft_threshold(sft, output_path)

  if (missing(Power)) stop("Choose a power based on softthresholdingpower.pdf")

  # Network construction
  net <- blockwiseModules(datExpr, power = Power, TOMType = "signed", minModuleSize = 30,
                          reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE,
                          pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "Salva", verbose = 3)

  # Save module info and plots
  rownames(datTraits) <- datTraits$Sample
  datTraits$Sample <- NULL
  save_module_info(net, datExpr, datTraits, output_path, Name)
  plot_module_trait_relationships(net, datExpr, datTraits, output_path)
  plot_module_summaries(net, TPM_filt, sample_table_some, Variables, Colors_plot, NumberCol, Name, output_path)

  # GO enrichment
  if (!dir.exists("topGO")) {
    dir.create("topGO", showWarnings = FALSE)
  }

  topGO_path <- file.path(output_path, "topGO")
  topGO_results <- run_topGO_for_modules(geneID2GO = geneID2GO, Name = Name, input_path = output_path, output_path_topGO = topGO_path, ontology = ontology,
                                         algorithm = algorithm, statistic = statistic, plot_similarity = plot_similarity,
                                         Number_GOs = Number_GOs,orgdb = orgdb, semdata = semdata, save_GeneNames = save_GeneNames, Ontologies = Ontologies, Annotation = Annotation)

  return(topGO_results)
}
