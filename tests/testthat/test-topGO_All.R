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
  geneUniverse <- paste0("Gene", 1:100)
  DEG <- sample(geneUniverse, 10)
  gene_factor <- prepare_topGO_genes(DEG, geneUniverse)
  geneID2GO <- setNames(replicate(100, sample(c("GO:0008150", "GO:0003674"), 2), simplify = FALSE), geneUniverse)
  result <- run_topGO_analysis(gene_factor, geneID2GO)
  expect_true("GOdata" %in% names(result))
  expect_true("result" %in% names(result))
})

test_that("process_topGO_results returns a data.frame or NULL", {
  geneUniverse <- paste0("Gene", 1:100)
  DEG <- sample(geneUniverse, 10)
  gene_factor <- prepare_topGO_genes(DEG, geneUniverse)
  geneID2GO <- setNames(replicate(100, sample(c("GO:0008150", "GO:0003674"), 2), simplify = FALSE), geneUniverse)
  analysis <- run_topGO_analysis(gene_factor, geneID2GO)
  results <- process_topGO_results(analysis$GOdata, analysis$result)
  expect_true(is.null(results) || is.data.frame(results))
})

test_that("topGO_All runs end-to-end", {
  geneUniverse <- paste0("Gene", 1:100)
  DEG <- sample(geneUniverse, 10)
  geneID2GO <- setNames(replicate(100, sample(c("GO:0008150", "GO:0003674"), 2), simplify = FALSE), geneUniverse)
  result <- topGO_All(DEG, geneID2GO, name = "test", output_dir = tempdir(), plot_similarity = FALSE)
  expect_true(is.list(result))
  expect_true("results_table" %in% names(result))
})
