#' Prepare genes for topGO analysis
#'
#' This function takes a gene vector or a data.frame with differential expression results
#' and returns a binary factor that shows which genes are in the list of interest.
#'
#' @param DEG Vector of genes or data.frame with differential expression results. Genes must be rownames in the data.frame.
#' @param geneUniverse Vector with all gene names of the organism universe (e.g., `geneID2GO` names).
#' @param padj_threshold Significance threshold for `padj` if `DEG` is a data.frame.
#'
#' @return Binary factor with gene names, to be used in `topGO`.
#' @export
prepare_topGO_genes <- function(DEG, geneUniverse, padj_threshold = 0.05) {
  if (is.vector(DEG)) {
    selected_genes <- DEG
  } else if (is.data.frame(DEG)) {
    if (!"padj" %in% colnames(DEG)) {
      warning("No 'padj' column found. All genes will be considered.")
      selected_genes <- rownames(DEG)
    } else {
      selected_genes <- rownames(DEG)[DEG$padj <= 0.05]
      if (is.null(selected_genes)) {
        stop("None of the selected genes are below the padj threshold.")
      }
    }
  } else {
    stop("None of the selected genes match the provided gene universe. Check gene IDs and annotation.")
  }

  if (!any(selected_genes %in% geneUniverse)) {
    stop("Gene IDs must correspond with the annotation file provided")
  }
  gene_factor <- factor(as.integer(geneUniverse %in% selected_genes))
  names(gene_factor) <- geneUniverse
  return(gene_factor)
}

#' Save topGO results to a text file
#'
#' @param go_results Data.frame with GO results.
#' @param output_path Path to save the file.
#' @param simplified Logical. If TRUE, only p-values are saved.
#' @export
save_topGO_results <- function(go_results, output_path, simplified = FALSE) {
  if (simplified) {
    utils::write.table(go_results %>% dplyr::select(pval),
                       file = output_path, sep = "\t", quote = FALSE,
                       row.names = TRUE, col.names = FALSE)
  } else {
    utils::write.table(go_results,
                       file = output_path, sep = "\t", quote = FALSE,
                       row.names = TRUE)
  }
}


#' Run GO enrichment analysis using topGO
#'
#' This function creates a topGOdata object and runs the enrichment test.
#'
#' @param gene_factor Binary factor indicating genes of interest (output of `prepare_topGO_genes()`).
#' @param geneID2GO Named list mapping gene IDs to GO terms.
#' @param ontology GO ontology to use: "BP", "MF", or "CC".
#' @param algorithm Algorithm for topGO test (default: "weight01").
#' @param statistic Statistical test to use (default: "fisher").
#'
#' @return A list with the topGOdata object and the result of the test.
#' @export
run_topGO_analysis <- function(gene_factor, geneID2GO,
                               ontology = "BP",
                               algorithm = "weight01",
                               statistic = "fisher") {
  if (!requireNamespace("topGO", quietly = TRUE)) {
    stop("The 'topGO' package is required but not installed.")
  }

  if (!is.factor(gene_factor)) {
    stop("gene_factor must be a factor.")
  }

  if (!is.list(geneID2GO)) {
    stop("geneID2GO must be a named list.")
  }

  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = gene_factor,
                annot = topGO::annFUN.gene2GO,
                gene2GO = geneID2GO)

  result <- topGO::runTest(GOdata, algorithm = algorithm, statistic = statistic)

  return(list(GOdata = GOdata, result = result))
}


#' Process topGO enrichment results
#'
#' This function filters significant GO terms, calculates enrichment, and returns a tidy result table.
#'
#' @param GOdata A `topGOdata` object.
#' @param result A result object from `runTest()`.
#' @param pval_threshold P-value cutoff to consider GO terms as significant.
#' @param gene_count Total number of genes tested (used for enrichment calculation).
#'
#' @return A data.frame with GO term statistics, p-values, and enrichment scores.
#' @export
process_topGO_results <- function(GOdata, result, pval_threshold = 0.05, gene_count = NULL) {
  if (!requireNamespace("topGO", quietly = TRUE)) {
    stop("The 'topGO' package is required but not installed.")
  }

  # Get significant GO terms
  pvals <- topGO::score(result)
  sig_terms <- pvals[pvals <= pval_threshold]

  if (length(sig_terms) == 0) {
    warning("No significant GO terms found.")
    return(NULL)
  }

  # Get term statistics
  stats <- topGO::termStat(GOdata, whichGO = names(sig_terms))
  stats$Term <- Term(rownames(stats))
  stats$pval <- sig_terms

  # If number of genes in list is not provided, get annotated genes from topGO object
  if (is.null(gene_count)) {
    gene_count <- length(topGO::sigGenes(GOdata))
  }

  total_genes <- length(unique(GOdata@allGenes))
  stats$Enrichment <- (stats$Significant / stats$Annotated) / (gene_count / total_genes)

  # Filter and order
  stats <- stats %>%
    dplyr::filter(Significant > 1, Enrichment > 1) %>%
    dplyr::arrange(Enrichment, Significant) %>%
    dplyr::mutate(Term = factor(Term, levels = Term))

  return(stats)
}


#' Plot topGO enrichment results
#'
#' This function creates a bubble plot of enriched GO terms using ggplot2.
#'
#' @param go_results A data.frame with GO enrichment results (output of `process_topGO_results()`).
#' @param title Title of the plot.
#' @param output_path Optional path to save the plot as PDF. If NULL, the plot is not saved.
#' @param units Units to set plot dimensions. Default to "in" (inches)
#' @param width Width units for the plot. Default to 10.
#' @param height Height units for the plot. Default to max(4, round(8 * nrow(go_results) / 50)), which was empirically obtained to get a nice heigh depending on the number of enriched GO terms.
#'
#' @return A ggplot object.
#' @export
plot_topGO_results <- function(go_results, title = "GO Enrichment", output_path = NULL, units = "in",
                               width = 10,
                               height = max(4, round(8 * nrow(go_results) / 50))) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The 'ggplot2' package is required but not installed.")
  }

  if (is.null(go_results) || nrow(go_results) == 0) {
    warning("No GO results to plot.")
    return(NULL)
  }

  go_results$Term <- factor(go_results$Term, levels = rev(go_results$Term))

  p <- ggplot2::ggplot(go_results, ggplot2::aes(x = Enrichment, y = Term, color = pval)) +
    ggplot2::geom_point(ggplot2::aes(size = Significant)) +
    ggplot2::scale_color_gradient2(low = "#253494", mid = "#41b6c4", high = "#edf8b1", midpoint = 0.025) +
    ggplot2::scale_size_binned() +
    ggplot2::labs(title = title, x = "Enrichment", y = NULL) +
    ggplot2::theme_gray() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "#d4d4d4", linetype = "dashed"),
      axis.text.y = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0, face = "bold")
    )

  if (!is.null(output_path)) {
    ggplot2::ggsave(output_path, plot = p, device = "pdf", width = width,
                    height = height,
                    units = units, scale = 1.5, limitsize = FALSE)
  }

  return(p)
}


#' Analyze GO term similarity and visualize clusters
#'
#' This function calculates semantic similarity between GO terms, reduces redundancy,
#' and generates scatterplot and treemap visualizations.
#'
#' @param go_results A data.frame with GO enrichment results (must include `pval`).
#' @param orgdb OrgDb package name (e.g., "org.At.tair.db").
#' @param ontology GO ontology used ("BP", "MF", or "CC"). Defaults to "BP".
#' @param semdata Precomputed semantic data (optional).
#' @param output_prefix Optional prefix for saving plots (PDFs).
#' @param Number_GOs Number of top GO term names to plot in the scatterplot. Defaults to 20.
#' @param units Units to set plot dimensions. Default to "in" (inches)
#' @param width Width units for the plot. Default to 10.
#' @param height Height units for the plot. Default to max(4, round(8 * nrow(go_results) / 50)), which was empirically obtained to get a nice heigh depending on the number of enriched GO terms.

#'
#' @return A list with the similarity matrix, reduced terms, and ggplot objects.
#' @export
analyze_GO_similarity <- function(go_results, orgdb = "org.At.tair.db",
                                  ontology = "BP", semdata = NULL,
                                  output_prefix = NULL, Number_GOs = 20,
                                  units = "in", width = 10, height = 8) {
  if (!requireNamespace("GOSemSim", quietly = TRUE) ||
      !requireNamespace("rrvgo", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("ggrepel", quietly = TRUE) ||
      !requireNamespace("treemap", quietly = TRUE) ||
      !requireNamespace("factoextra", quietly = TRUE)) {
    stop("Required packages: GOSemSim, rrvgo, ggplot2, ggrepel, treemap, factoextra.")
  }

  terms <- rownames(go_results)
  scores <- setNames(-log10(go_results$pval), terms)

  if (is.null(semdata)) {
    semdata <- GOSemSim::godata(orgdb, ont = ontology)
  }

  simMatrix <- rrvgo::calculateSimMatrix(terms, orgdb = orgdb,
                                         ont = ontology, semdata = semdata,
                                         method = "Rel")

  if (nrow(simMatrix) < 2) {
    warning("Not enough terms for similarity analysis.")
    return(NULL)
  }

  reducedTerms <- rrvgo::reduceSimMatrix(simMatrix, scores = scores,
                                         threshold = 0.5, orgdb = orgdb)

  coords <- stats::cmdscale(as.dist(1 - simMatrix), k = 2, eig=TRUE)$points
  df <- cbind(as.data.frame(coords),
              reducedTerms[match(rownames(coords), reducedTerms$go),
                           c("term", "parent", "parentTerm", "score")]) %>%
    dplyr::arrange(-score)

  top_terms <- df$term[1:min(Number_GOs, nrow(df))]

  # Clustering
  kmax <- if (nrow(df) <= 10) nrow(df) - 1 else 10
  viz <- factoextra::fviz_nbclust(df[, 1:2], kmeans, method = "silhouette", k.max = kmax)
  nclust <- as.integer(viz$data[viz$data$y == max(viz$data$y),"clusters"])
  set.seed(123)
  df$clust <- as.factor(kmeans(df[, 1:2], centers = nclust)$cluster)
  df$label <- ifelse(df$term %in% top_terms, df$parentTerm, "")

  df_modified <- df %>%
    tibble::rownames_to_column("GO_ID") %>%
    dplyr::group_by(clust) %>%
    dplyr::group_modify(~ {
      # .x contains current clúster data
      data_current_cluster <- .x

      # Verify if all labels in this cluster are empty
      if (all(data_current_cluster$label == "")) {
        # If all labels are empty:
        # 1. Arrange in desc according to score values
        cluster_modified <- data_current_cluster %>%
          arrange(desc(score))

        # 2. Identify rows to be updated (first two, unless cluster is smaller)
        rows_to_update_idx <- head(1:nrow(cluster_modified), 2)

        # 3. Update 'label' column of these rows with 'partenTerm' value
        if (length(rows_to_update_idx) > 0) {
          cluster_modified$label[rows_to_update_idx] <- cluster_modified$parentTerm[rows_to_update_idx]
        }
        # Return modified cluster (and arranged according to 'score')
        return(cluster_modified)
      } else {
        # Return previous cluster if there were labelled rows already
        return(data_current_cluster)
      }
    }) %>%
    ungroup() %>% # Ungroup to obtain final modified data frame
    tibble::column_to_rownames("GO_ID")


  # Scatterplot
  scatterplot <- ggplot2::ggplot(df_modified, ggplot2::aes(x = V1, y = V2, color = clust)) +
    ggplot2::geom_point(ggplot2::aes(size = score), alpha = 0.5) +
    ggplot2::scale_color_manual(guide = "none", values = as.vector(ggsci::pal_npg("nrc")(10))) +
    ggplot2::scale_size_continuous(guide = "none", range = c(0, 25)) +
    ggplot2::scale_x_continuous(name="") +
    ggplot2::scale_y_continuous(name="") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                   axis.text.y=ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())+
    ggrepel::geom_label_repel(ggplot2::aes(label = label),
                              data = subset(df_modified, parent == rownames(df_modified)),
                              box.padding = grid::unit(1, "lines"),
                              size = 10, max.overlaps = 100)



  # Treemap
  treemap_file <- NULL
  if (!is.null(output_prefix)) {
    scatter_file <- paste0(output_prefix, "_Scatterplot.pdf")
    treemap_file <- paste0(output_prefix, "_Treemap.pdf")
    ggplot2::ggsave(scatter_file, plot = scatterplot, width = 18, height = 10, units = units, scale= 1.5, limitsize = FALSE)

    grDevices::pdf(treemap_file, width = width, height = height)
    treemap::treemap(reducedTerms,
                     index = c("parentTerm", "term"),
                     vSize = "score",
                     type = "index",
                     title = "",
                     palette = rep(ggsci::pal_npg("nrc")(10), length(unique(reducedTerms$parent))),
                     fontcolor.labels = c("#FFFFFFDD", "#00000080"),
                     bg.labels = 0,
                     border.col = "#00000080")
    grDevices::dev.off()
  }

  return(list(similarity = simMatrix,
              reduced = reducedTerms,
              scatterplot = scatterplot,
              treemap_file = treemap_file))
}


#' GeneNames_GOs: Generate Excel report of genes associated with GO terms
#'
#' This function generates an Excel report of genes associated with specified GO terms.
#' It reads a functional annotation file, filters genes based on the provided DEGs and GO terms,
#' and saves the results in an Excel file.
#'
#' @param Annotation Two column data frame. First column must include Gene IDs.
#' Second column must include functional annotation.
#' @param geneID2GO Named list mapping gene IDs to GO terms.
#' @param output_dir Directory where the output file will be saved.
#' @param DEG Vector or Data frame. If vector, it must include Gene IDs of interest.
#' If Data frame, it must be the output of differential expressed genes (row names should be Gene IDs).
#' @param Ontologies Character vector of GO terms to search for.
#' @param name Base name for the output file.
#'
#' @return No return value. An Excel file is written to disk.
#' @export
GeneNames_GOs <- function(Annotation, geneID2GO, output_dir, DEG, Ontologies = NULL, name) {
  if (is.vector(DEG)) {
    genes <- DEG
  } else if (is.data.frame(DEG)) {
    genes <- rownames(DEG)
  } else {
    stop("You should provide a vector or a data frame as `DEG` argument.")
  }

  geneID2GO_stack <- stack(geneID2GO)
  if (!any(genes %in% geneID2GO_stack$ind)) {
    warning("Gene IDs do not match with the annotation geneID2GO. The output will be empty.")
  }

  geneID2GO_stack_filtered <- geneID2GO_stack[geneID2GO_stack$ind %in% genes,]

  if (!is.null(Ontologies)) {
    if (any(Ontologies %in% geneID2GO_stack_filtered$values)) {
      geneID2GO_stack_filtered <- geneID2GO_stack_filtered[geneID2GO_stack_filtered$values %in% Ontologies,]
    } else {
      stop("Any of your GO terms of interest is representing Gene IDs of this analysis.")
    }
  }
  List_GOs_Genes <- split(geneID2GO_stack_filtered$ind, geneID2GO_stack_filtered$values)

  if (ncol(Annotation) < 2) {
    stop("`Annotation` must have at least two columns: Gene ID and Description.")
  }

  if (!any(genes %in% Annotation[,1])) {
    warning("Gene IDs do not match with the `Annotation` argument. No description will be added.")
  }

  # Create final report
  Final_Report <- list()
  for (GO in names(List_GOs_Genes)) {
    gene_df <- data.frame(Gene = List_GOs_Genes[[GO]]) %>%
      dplyr::left_join(Annotation)
    Term <- stringr::str_sub(AnnotationDbi::Term(GO), 1, 31)
    Final_Report[[Term]] <- gene_df
  }

  # Write to Excel
  openxlsx::write.xlsx(Final_Report,
                       file = file.path(output_dir, paste0(name, "_Ontologies_Genes.xlsx")))
  message("Excel file saved to: ", file.path(output_dir, paste0(name, "_Ontologies_Genes.xlsx")))
  return(Final_Report)
}




#' Run full GO enrichment analysis with topGO
#'
#' This function performs the full GO enrichment workflow: gene preparation, enrichment test,
#' result processing, visualization, and optional ancestor/similarity analysis.
#'
#' @param DEG Vector or data.frame of differential expressed genes.
#' @param geneID2GO Named list mapping gene IDs to GO terms.
#' @param name Identifier for output files. If genes comes from a biological comparison, it is recommended to use that name
#' @param output_dir Directory to save results (default: current working directory).
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
#' @return A list with all intermediate and final results.
#' @export
topGO_All <- function(DEG, geneID2GO, name = "GO_analysis",
                      output_dir = ".", ontology = "BP",
                      algorithm = "weight01",
                      statistic = "fisher",
                      plot_similarity = TRUE, Number_GOs = 20,
                      orgdb = "org.At.tair.db", semdata = NULL,
                      save_GeneNames = FALSE, Annotation, Ontologies = NULL) {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  output_prefix <- file.path(output_dir, paste0("topGO_", name))

  # 1. Prepare genes
  geneUniverse <- names(geneID2GO)
  gene_factor <- prepare_topGO_genes(DEG, geneUniverse)

  # 2. Execute GO analysis
  go_analysis <- run_topGO_analysis(gene_factor, geneID2GO, ontology = ontology, algorithm = algorithm, statistic = statistic)

  # 3. Process results
  go_results <- process_topGO_results(go_analysis$GOdata, go_analysis$result,
                                      gene_count = sum(gene_factor == 1))

  if (is.null(go_results)) {
    warning("No significant GO terms found.")
    return(NULL)
  }

  # 4. Save results
  save_topGO_results(go_results, paste0(output_prefix, "_results.txt"))
  save_topGO_results(go_results, paste0(output_prefix, "_pvals.txt"), simplified = TRUE)

  # 5. Principal plot
  bubble_plot <- plot_topGO_results(go_results, title = name,
                                    output_path = paste0(output_prefix, "_bubble.pdf"))

  output_list <- list(
    GOdata = go_analysis$GOdata,
    result = go_analysis$result,
    results_table = go_results,
    bubble_plot = bubble_plot
    )

  # 6. Semantic plots (scatterplot and treemap)
  if (plot_similarity) {
    similarity_results <- analyze_GO_similarity(go_results, orgdb = orgdb,
                                                ontology = ontology,
                                                semdata = semdata, Number_GOs = Number_GOs,
                                                output_prefix = output_prefix)
    output_list$similarity_results <- similarity_results
  }

  if (save_GeneNames){
    GeneNames <- GeneNames_GOs(Annotation = Annotation, geneID2GO = geneID2GO,
                               output_dir = output_dir, DEG = DEG,
                               Ontologies = Ontologies, name = name)
    output_list$GeneNames <- GeneNames
  }

  return(output_list)
}

