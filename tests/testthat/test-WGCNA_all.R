
library(testthat)
library(RNAseqEasy)
library(WGCNA)
library(topGO)

test_that("WGCNA_Modules runs without errors and generates expected output files", {
  # Simulate sample_table
  sample_table <- data.frame(
    Sample = paste0("Wendy", seq(1,12)),
    Genotype = rep(c("Tak1", "gh3a"), each = 6),
    Treatment = rep(c("Mock", "OPDA"), each = 3),
    Replicate = seq(1,3)
  )

  data("output_DESeq2")

  data("Marchantia7_tx2gene")

  data("Mpo_GO_GOSLIM")
  geneID2GO <- load_topGO_db(Mpo_GO_GOSLIM, "Transcript", "GO")

  sampleDir <- "/Users/salva/Google Drive/My Drive/Antiguo Drive/CNB/Ayuda Wendy/RNAseq/02_Salmon/"

  output_path <- "/Users/salva/Desktop/test_RNAseqEasy/"

  save <- GOSemSim::godata("org.At.tair.db", ont="BP")


  # Run WGCNA_Modules
  result <- WGCNA_Modules(
    output_path = output_path,
    sampleDir = sampleDir,
    sample_table = sample_table,
    DEGs = rownames(output_DESeq2),
    Variables = c("Genotype", "Treatment"),
    tx2gene = Marchantia7_tx2gene,
    Power = 18,
    Name = "Test",
    Colors_plot = NULL,
    NumberCol = 1,
    geneID2GO = geneID2GO,
    semdata = save
  )

  # Check that the expected output files are generated
  expect_true(file.exists(file.path(output_path, "Test_geneInfo.csv")))
  expect_true(file.exists(file.path(output_path, "sampleClustering.pdf")))
  expect_true(file.exists(file.path(output_path, "Summary_Modules_Test.pdf")))
})
