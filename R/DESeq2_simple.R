#' Save DESeq2 results to text files
#'
#' This function writes the full and significant DESeq2 results to tab-delimited text files.
#'
#' @param res A DESeqResults object.
#' @param output_dir Directory where the results will be saved.
#' @param name Base name for the output files.
#'
#' @return No return value. Files are written to disk.
#' @export
save_deseq_results <- function(res, output_dir, name) {
  write.table(res, file = file.path(output_dir, paste0(name, ".txt")),
              quote = FALSE, sep = "	")
  write.table(subset(res, padj <= 0.05),
              file = file.path(output_dir, paste0(name, "_Sig.txt")),
              quote = FALSE, sep = "	")
}



#' Plot heatmap of significant DEGs
#'
#' This function generates a heatmap of genes with adjusted p-value <= 0.05 using normalized counts.
#'
#' @param dds A DESeqDataSet object.
#' @param res A DESeqResults object.
#' @param variables Character vector of column names in `sample_table` to annotate samples in heatmap.
#' @param name Name used for saving the plot.
#' @param width Width inches for the plot. Default to 6.
#' @param height Height inches for the plot. Default to 8.
#' @param output_dir Directory where the plot will be saved.
#'
#' @return No return value. A PDF heatmap is saved.
#' @export
plot_deseq_heatmap <- function(dds, res, variables, name, output_dir, width = 6, height = 8) {
  f1 <- function(x) stringr::str_split(x, pattern = "_")[[1]][1] # In each sample name there are different things separated by "_" (name, run, etc...). With this function we will get the first one, that is, the real sample name

  select <-  as.data.frame(subset(res, padj <= 0.05)) # we select the genes we want to plot
  ntd <- normTransform(dds) # this gives log2(n + 1)
  heatmap_data <- assay(ntd)[rownames(select),]
  colnames(heatmap_data) <- sapply(colnames(heatmap_data), f1)
  heatmap_data_ann <- as.data.frame(colData(dds)) %>% dplyr::select(matches(variables)) # we create a data-frame with the phenotypic labels to plot
  rownames(heatmap_data_ann) <- sapply(rownames(heatmap_data_ann), f1)

  # Create heatmap
  pheatmap::pheatmap(heatmap_data, annotation_col = heatmap_data_ann,
                     show_rownames = FALSE, fontsize_row = 8, fontsize_col = 8,
                     filename = file.path(output_dir, paste0(name, "_heatmap.pdf")),
                     width = width, height = height, cluster_rows=TRUE, cluster_cols=TRUE,
                     scale = "row")
}



#' Run DESeq2 analysis with optional LRT test and GO enrichment
#'
#' This function performs differential expression analysis using DESeq2, optionally
#' including a likelihood ratio test (LRT) with a reduced design. It also generates
#' heatmaps and runs GO enrichment analysis for significant genes.
#'
#' @param output_path Directory where results will be saved.
#' @param sampleDir Directory containing sample folders with quant.sf files.
#' @param sample_table Data frame with sample metadata.
#' @param Include A vector of condition(s) to include samples from `sample_table`.
#' It will include samples matching this/these condition(s) and exclude all others.
#' Default to NULL.
#' @param Exclude A vector of condition(s) to exclude samples from `sample_table`.
#' Default to NULL.
#' @param Variable Character. Column name in sample_table used for grouping.
#' @param Design Character. Design formula for DESeq2 (e.g., "~ condition").
#' @param Group Character. Whether to use custom contrasts ("YES" or "NO").
#' @param Name Character. Base name for output files.
#' @param Contrast Character Custom contrasts for DESeq2 results.
#' @param Reduced Logical. Whether to use a reduced design (LRT test). Default is FALSE.
#' @param Reduced_design A formula specifying the reduced model (required if `reduced = TRUE`).
#' @param log2FCtopGO Numeric. log2(fold change) threshold for GO analysis. Default is 1.
#' @param geneID2GO Named list mapping gene IDs to GO terms.
#' @param ontology GO ontology to use ("BP", "MF", "CC"). Defaults to "BP".
#' @param algorithm Algorithm for topGO test (default: "weight01").
#' @param statistic Statistical test to use (default: "fisher").
#' @param plot_similarity Logical, whether to analyze and visualize GO term similarity.
#' @param Number_GOs Number of top GO term names to plot in the scatterplot. Defaults to 20.
#' @param orgdb OrgDb package name for similarity analysis. Defaults to "org.At.tair.db".
#' @param semdata Optional precomputed semantic data.
#'
#' @return No return value. Results are saved to disk.
#' @export
DESeq2_simple <- function(output_path, sampleDir, sample_table, Include = NULL, Exclude = NULL,
                          tx2gene, min_count = 10, min_samples = 3,
                          Variable, Design, Group = "NO", Name, Contrast,
                          Reduced = FALSE, Reduced_design = NULL,
                          log2FCtopGO = 1, ontology = "BP", plot_similarity = TRUE,
                          geneID2GO, algorithm = "weight01",
                          statistic = "fisher", Number_GOs = 20,
                          orgdb = "org.At.tair.db", semdata = NULL) {

  sample_table_some <- get_sample_subset(sample_table, Include = Include,
                                         Exclude = Exclude)
  sample_table_some <- add_sample_path(sampleDir = sampleDir, sample_table = sample_table_some)
  txi <- load_tximport_data(sampleDir, sample_table_some, tx2gene)

  dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = sample_table_some, design = as.formula(paste("~", Design)))
  dds <- filter_low_counts(dds, min_count = min_count, min_samples = min_samples)

  if (Reduced) {
    dds <- DESeq2::DESeq(dds, test = "LRT", reduced = Reduced_design)
  } else {
    dds <- DESeq2::DESeq(dds)
  }

  if (Group == "YES"){
    First <- Variable[1]
    Second <- Variable[2]
    Levels <- as.character()
    for (first in levels(sample_table_some[,First])[levels(sample_table_some[,First]) %in% unique(sample_table_some[, First])]) {
      for (second in levels(sample_table_some[,Second])[levels(sample_table_some[,Second]) %in% unique(sample_table_some[, Second])]) {
        Levels <- append(Levels, paste0(first, second))
      }
    }
    dds$group <- factor(paste0(sample_table_some[,First], sample_table_some[,Second]), levels = Levels)
    design(dds) <- ~group
    dds <- DESeq(dds)
  }

  # get the model matrix
  mod_mat <- model.matrix(design(dds), colData(dds))

  Alt_contrasts <- as.character()

  # Define coefficient vectors for each condition
  for (Condition in Variable) {
    for (Element in unique(colData(dds)[,Condition])) {
      assign(Element, colMeans(mod_mat[colData(dds)[,Condition] == Element, ]), envir = .GlobalEnv)
      Alt_contrasts <- append(Alt_contrasts, Element)
    }
  }

  if (length(Variable) == 2) {
    for (Primero in unique(colData(dds)[,Variable[1]])) {
      for (Segundo in unique(colData(dds)[,Variable[2]])) {
        assign(paste(Primero, Segundo, sep = "_"), colMeans(mod_mat[colData(dds)[,Variable[1]] == Primero & colData(dds)[,Variable[2]] == Segundo, ]), envir = .GlobalEnv)
        Alt_contrasts <- append(Alt_contrasts, paste(Primero, Segundo, sep = "_"))
      }
    }
  }

  if (length(Variable) == 3) {
    for (Primero in unique(colData(dds)[,Variable[1]])) {
      for (Segundo in unique(colData(dds)[,Variable[2]])) {
        assign(paste(Primero, Segundo, sep = "_"), colMeans(mod_mat[colData(dds)[,Variable[1]] == Primero & colData(dds)[,Variable[2]] == Segundo, ]), envir = .GlobalEnv)
        Alt_contrasts <- append(Alt_contrasts, paste(Primero, Segundo, sep = "_"))
        for (Tercero in unique(colData(dds)[,Variable[3]])) {
          assign(paste(Primero, Segundo, Tercero, sep = "_"), colMeans(mod_mat[colData(dds)[,Variable[1]] == Primero & colData(dds)[,Variable[2]] == Segundo & colData(dds)[,Variable[3]] == Tercero, ]), envir = .GlobalEnv)
          Alt_contrasts <- append(Alt_contrasts, paste(Primero, Segundo, Tercero, sep = "_"))
        }
      }
      for (Tercero in unique(colData(dds)[,Variable[3]])) {
        assign(paste(Primero, Tercero, sep = "_"), colMeans(mod_mat[colData(dds)[,Variable[1]] == Primero & colData(dds)[,Variable[3]] == Tercero, ]), envir = .GlobalEnv)
        Alt_contrasts <- append(Alt_contrasts, paste(Primero, Tercero, sep = "_"))
      }
    }
  }


  if (missing(Contrast)) {
    print("These are the DESeq2 contrasts, that you have to use like list(c(\"ContrastA\", \"ContrastB\")) : ", quote= FALSE)
    print(resultsNames(dds))
    print("These are the costumized contrasts, that you have to use like ContrastA - ContrastB : ", quote= FALSE)
    print(Alt_contrasts)
    stop("Choose one of the previous contrasts")
  }


  # Guardar resultados
  res <- DESeq2::results(dds, contrast = Contrast)
  save_deseq_results(res, output_path, Name)

  # Heatmap
  plot_deseq_heatmap(dds, res, Variable, Name, output_path)

  # topGO
  up_genes <- rownames(subset(res, padj <= 0.05 & log2FoldChange >= log2FCtopGO))
  down_genes <- rownames(subset(res, padj <= 0.05 & log2FoldChange <= -log2FCtopGO))
  up_down_genes <- c(up_genes, down_genes)

  if (!dir.exists("topGO")) {
    dir.create("topGO", showWarnings = FALSE)
  }

  topGO_path <- file.path(output_path, "topGO")

  if (length(up_down_genes) > 1) {
    topGO_All(up_down_genes, geneID2GO = geneID2GO, name = paste0(Name, "_All"),
              output_dir = topGO_path, ontology = ontology,
              plot_similarity = plot_similarity, algorithm = algorithm, statistic = statistic,
              Number_GOs = Number_GOs, orgdb = orgdb, semdata = semdata)
  }

  if (length(up_genes) > 1) {
    topGO_Up <- topGO_All(up_genes, geneID2GO = geneID2GO, name = paste0(Name, "_Up"),
                          output_dir = topGO_path, ontology = ontology,
                          plot_similarity = plot_similarity, algorithm = algorithm, statistic = statistic,
                          Number_GOs = Number_GOs, orgdb = orgdb, semdata = semdata)
  }

  if (length(down_genes) > 1) {
    topGO_Down <- topGO_All(down_genes, geneID2GO = geneID2GO, name = paste0(Name, "_Down"),
                            output_dir = topGO_path, ontology = ontology,
                            plot_similarity = plot_similarity, algorithm = algorithm, statistic = statistic,
                            Number_GOs = Number_GOs, orgdb = orgdb, semdata = semdata)
  }

  Output <- list(Up = topGO_Up,
                 Down = topGO_Down)
  return(Output)
}

