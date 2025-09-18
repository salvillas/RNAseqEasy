library(testthat)
library(dplyr)
library(tximport)

# Dummy data for testing
sample_table <- data.frame(
  Folder = c("sample1", "sample2", "sample3"),
  Condition = c("A", "B", "A"),
  stringsAsFactors = FALSE
)

selection_yes <- "sample1"
selection_no <- "sample2"

test_that("get_sample_subset returns full table when Include of Exclude are NULL", {
  expect_equal(get_sample_subset(sample_table), sample_table)
  expect_equal(get_sample_subset(sample_table, Include = NULL, Exclude = NULL), sample_table)
})

test_that("get_sample_subset filters correctly with Include", {
  result <- get_sample_subset(sample_table, Include = selection_yes)
  expect_true(all(grepl("sample1", result$Folder)))
})

test_that("get_sample_subset filters correctly with Exclude", {
  result <- get_sample_subset(sample_table, Exclude = selection_no)
  expect_false(any(grepl("sample2", result$Folder)))
})


test_that("add_sample_path correctly adds Folder column", {

  sampleDir <- "/Users/salva/Google Drive/My Drive/Antiguo Drive/CNB/Ayuda Wendy/RNAseq/02_Salmon/"
  sample_table <- read.delim("/Users/salva/Google Drive/My Drive/Antiguo Drive/CNB/Ayuda Wendy/RNAseq/Sample_Data_Wendy.txt")

  # Run the function
  result <- add_sample_path(sampleDir, sample_table)

  # Check that Folder column exists
  expect_true("Folder" %in% colnames(result))

  # Check that all paths are correctly assigned
  expected_paths <- list.dirs(sampleDir, full.names = FALSE, recursive = FALSE)
  expect_equal(sort(result$Folder), sort(expected_paths))

})


test_that("add_sample_path throws error if no matching column is found", {
  sampleDir <- "/Users/salva/Google Drive/My Drive/Antiguo Drive/CNB/"

  sample_table <- data.frame(
    ID = c("X", "Y", "Z"),
    Group = c("A", "B", "A"),
    stringsAsFactors = FALSE
  )

  expect_error(
    add_sample_path(sampleDir, sample_table),
    "No column in 'sample_table' matches folder names"
  )
})


test_that("load_tximport_data returns a list with expected structure", {
  data("Marchantia7_tx2gene")

  sampleDir <- "/Users/salva/Google Drive/My Drive/Antiguo Drive/CNB/Ayuda Wendy/RNAseq/02_Salmon/"
  sample_table <- read.delim("/Users/salva/Google Drive/My Drive/Antiguo Drive/CNB/Ayuda Wendy/RNAseq/Sample_Data_Wendy.txt")

  load_tximport_data(sampleDir, sample_table, Marchantia7_tx2gene)
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
