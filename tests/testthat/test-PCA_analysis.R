library(testthat)
library(DESeq2)
library(ggplot2)

# Create a dummy DESeqDataSet for testing
set.seed(123)
counts <- matrix(rnbinom(10000, mu = 1000, size = 1), ncol = 10)
colnames(counts) <- paste0("sample", 1:10)
rownames(counts) <- paste0("gene", 1:1000)

col_data <- data.frame(
  condition = rep(c("A", "B"), each = 5),
  batch = rep(c("X", "Y"), times = 5),
  row.names = colnames(counts)
)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = ~ condition)

test_that("calculate_pca returns expected structure", {
  result <- calculate_pca(dds, ntop = 50)
  expect_type(result, "list")
  expect_true(all(c("pca", "percentVar", "vsd") %in% names(result)))
  expect_s3_class(result$pca, "prcomp")
  expect_true(is.numeric(result$percentVar))
  expect_s4_class(result$vsd, "DESeqTransform")
})

test_that("plot_pca returns a ggplot object", {
  result <- calculate_pca(dds, ntop = 50)
  colors <- c(A = "#1f77b4", B = "#ff7f0e")
  p <- plot_pca(result, variables = c("batch", "condition"), colors = colors, components = c(1, 2), name = "test", output_dir = tempdir(), width = 10, height = 8, units = "in")
  expect_s3_class(p, "ggplot")
})

test_that("run_pca_analysis returns a ggplot object", {
  colors <- c(A = "#1f77b4", B = "#ff7f0e")
  p <- run_pca_analysis(dds, variables = c("batch", "condition"), colors = colors, components = c(1, 2), name = "test", output_dir = tempdir(), ntop = 50, width = 8, height = 5, units = "in")
  expect_s3_class(p, "ggplot")
})
