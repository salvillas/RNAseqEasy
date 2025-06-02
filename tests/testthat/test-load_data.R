library(testthat)
library(dplyr)
library(tximport)

# Dummy data for testing
sample_table <- data.frame(
  Folder = c("sample1", "sample2", "sample3"),
  Condition = c("A", "B", "A"),
  stringsAsFactors = FALSE
)

selection_yes <- list(YES = "sample1")
selection_no <- list(NO = "sample2")

test_that("get_sample_subset returns full table when selection is NULL or 'None'", {
  expect_equal(get_sample_subset(sample_table, NULL), sample_table)
  expect_equal(get_sample_subset(sample_table, "None"), sample_table)
})

test_that("get_sample_subset filters correctly with YES", {
  result <- get_sample_subset(sample_table, selection_yes)
  expect_true(all(grepl("sample1", result$Folder)))
})

test_that("get_sample_subset filters correctly with NO", {
  result <- get_sample_subset(sample_table, selection_no)
  expect_false(any(grepl("sample2", result$Folder)))
})

test_that("load_tximport_data returns a list with expected structure", {
  skip("Requires actual quant.sf files and tx2gene mapping to test properly.")
})

test_that("filter_low_counts filters genes correctly", {
  # Create a dummy DESeqDataSet
  counts_matrix <- matrix(c(10, 5, 0, 20, 30, 5, 0, 0, 0), nrow = 3)
  colnames(counts_matrix) <- c("sample1", "sample2", "sample3")
  rownames(counts_matrix) <- c("gene1", "gene2", "gene3")

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts_matrix,
    colData = data.frame(row.names = colnames(counts_matrix), condition = c("A", "B", "A")),
    design = ~ condition
  )

  filtered_dds <- filter_low_counts(dds, min_count = 10, min_samples = 2)
  expect_equal(nrow(filtered_dds), 1)
  expect_equal(rownames(filtered_dds), "gene1")
})
