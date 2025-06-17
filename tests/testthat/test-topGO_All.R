library(testthat)
library(topGO)
library(RNAseqEasy)

test_that("prepare_topGO_genes works with vector input", {
  geneUniverse <- paste0("Gene", 1:100)
  DEG <- sample(geneUniverse, 10)
  gene_factor <- prepare_topGO_genes(DEG, geneUniverse)
  expect_true(is.factor(gene_factor))
  expect_equal(length(gene_factor), length(geneUniverse))
  expect_equal(sum(gene_factor == 1), length(DEG))
})

test_that("prepare_topGO_genes works with data.frame input", {
  geneUniverse <- paste0("Gene", 1:100)
  df <- data.frame(padj = runif(100, max = 0.1), row.names = geneUniverse)
  gene_factor <- prepare_topGO_genes(df, geneUniverse, padj_threshold = 0.05)
  expect_true(is.factor(gene_factor))
})

test_that("run_topGO_analysis returns expected structure", {
  data("DEG_example_list")
  data("Mpo_GO_GOSLIM")
  geneID2GO <- load_topGO_db(Mpo_GO_GOSLIM, "Transcript", "GO")
  geneUniverse <- unique(Mpo_GO_GOSLIM$Transcript)
  gene_factor <- prepare_topGO_genes(DEG_example_list, geneUniverse)
  result <- run_topGO_analysis(gene_factor, geneID2GO)
  expect_true("GOdata" %in% names(result))
  expect_true("result" %in% names(result))
})

test_that("process_topGO_results returns a data.frame or NULL", {
  data("DEG_example_list")
  data("Mpo_GO_GOSLIM")
  geneID2GO <- load_topGO_db(Mpo_GO_GOSLIM, "Transcript", "GO")
  geneUniverse <- unique(Mpo_GO_GOSLIM$Transcript)
  gene_factor <- prepare_topGO_genes(DEG_example_list, geneUniverse)
  analysis <- run_topGO_analysis(gene_factor, geneID2GO)
  results <- process_topGO_results(analysis$GOdata, analysis$result)
  expect_true(is.null(results) || is.data.frame(results))
})

test_that("topGO_All runs end-to-end", {
  data("DEG_example_list")
  data("Mpo_GO_GOSLIM")
  geneID2GO <- load_topGO_db(Mpo_GO_GOSLIM, "Transcript", "GO")
  geneUniverse <- unique(Mpo_GO_GOSLIM$Transcript)
  gene_factor <- prepare_topGO_genes(DEG_example_list, geneUniverse)
  result <- topGO_All(DEG_example_list, geneID2GO, name = "test", output_dir = tempdir(), plot_similarity = FALSE)
  expect_true(is.list(result))
  expect_true("results_table" %in% names(result))
})

test_that("topGO_All runs end-to-end and saves GO terms of interest", {
  data("DEG_example_list")
  data("Mpo_GO_GOSLIM")
  data("Annotation_v6.1_Genes")
  geneID2GO <- load_topGO_db(Mpo_GO_GOSLIM, "Transcript", "GO")
  geneUniverse <- unique(Mpo_GO_GOSLIM$Transcript)
  gene_factor <- prepare_topGO_genes(DEG_example_list, geneUniverse)
  GOterms_of_interest <- c("GO:0009617", "GO:0042742", "GO:0002215")
  result <- topGO_All(DEG_example_list, geneID2GO, name = "test", output_dir = tempdir(), plot_similarity = FALSE,
                      save_GeneNames = TRUE, Annotation = Annotation_v6.1_Genes, Ontologies = GOterms_of_interest)
  expect_true(is.list(result))
  expect_true("results_table" %in% names(result))
})
