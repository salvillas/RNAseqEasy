library(testthat)
library(RNAseqEasy)
library(DESeq2)

test_that("TPM_all returns a data.frame, output file is created and data.frame is not empty", {
  sample_table <- data.frame(
    Sample = paste0("Wendy", seq(1,12)),
    Genotype = rep(c("Tak1", "gh3a"), each = 6),
    Treatment = rep(c("Mock", "OPDA"), each = 3),
    Replicate = seq(1,3)
  )
  data("Marchantia7_tx2gene")
  sampleDir <- "/Users/salva/Google Drive/My Drive/Antiguo Drive/CNB/Ayuda Wendy/RNAseq/02_Salmon/"
  output_path <- "/Users/salva/Desktop/test_RNAseqEasy/"
  TPM_Wendy <- TPM_all(sampleDir = sampleDir, sample_table = sample_table,
                       output_path = output_path, Include = NULL, Exclude = NULL,
                       Save_results = TRUE, tx2gene = Marchantia7_tx2gene)
  expect_true(is.data.frame(TPM_Wendy))
  expect_true(file.exists(file.path(output_path, "TPM_all_samples.txt")))
  expect_gt(nrow(TPM_Wendy), 0)
  expect_gt(ncol(TPM_Wendy), 0)
})


