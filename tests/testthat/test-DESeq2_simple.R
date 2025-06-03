
library(testthat)
library(RNAseqEasy)
library(topGO)

# Crear datos simulados mínimos
sample_table <- data.frame(
  Sample = paste0("Wendy", seq(1,12)),
  Genotype = rep(c("Tak1", "gh3a"), each = 6),
  Treatment = rep(c("Mock", "OPDA"), each = 3),
  Replicate = seq(1,3)
)
sample_table$Treatment <- factor(sample_table$Treatment, levels = c("Mock", "OPDA"))

data("Marchantia7_tx2gene")

data("Mpo_GO_GOSLIM")
geneID2GO <- load_topGO_db(Mpo_GO_GOSLIM, "Transcript", "GO")

sampleDir <- "/Users/salva/Google Drive/My Drive/Antiguo Drive/CNB/Ayuda Wendy/RNAseq/02_Salmon/"

output_path <- "/Users/salva/Desktop/test_RNAseqEasy/"

save <- GOSemSim::godata("org.At.tair.db", ont="BP")

test_that("DESeq2_simple stops when no Contrast is provided", {
  expect_error(result <- DESeq2_simple(
    output_path = output_path,
    sampleDir = sampleDir,
    sample_table = sample_table,
    Include = NULL,
    Exclude = NULL,
    tx2gene = Marchantia7_tx2gene,
    Variable = c("Genotype", "Treatment"),
    Design = "Genotype + Treatment + Genotype:Treatment",
    Group = "NO",
    Name = "test_DESeq2",
    Reduced = FALSE,
    log2FCtopGO = 1,
    geneID2GO = geneID2GO,
    ontology = "BP",
    plot_similarity = TRUE,
    orgdb = "org.At.tair.db",
    semdata = save
  ),
  "Choose one of the previous contrasts")
})


test_that("DESeq2_simple runs without errors and generates expected output files", {

  result <- DESeq2_simple(
    output_path = output_path,
    sampleDir = sampleDir,
    sample_table = sample_table,
    Include = NULL,
    Exclude = NULL,
    tx2gene = Marchantia7_tx2gene,
    Variable = c("Genotype", "Treatment"),
    Design = "Genotype + Treatment + Genotype:Treatment",
    Group = "NO",
    Name = "test_DESeq2",
    Contrast = list(c("Treatment_OPDA_vs_Mock")),
    Reduced = FALSE,
    log2FCtopGO = 1,
    geneID2GO = geneID2GO,
    ontology = "BP",
    plot_similarity = TRUE,
    orgdb = "org.At.tair.db",
    semdata = save
  )

  # Verificar que se generaron los archivos esperados
  expect_true(file.exists(file.path(output_path, "test_DESeq2.txt")))
  expect_true(file.exists(file.path(output_path, "test_DESeq2_Sig.txt")))
  expect_true(file.exists(file.path(output_path, "test_DESeq2_heatmap.pdf")))
})
