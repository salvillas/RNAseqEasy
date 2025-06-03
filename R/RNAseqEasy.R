library(usethis)
library(roxygen2)
library(devtools)
library(testthat)
library(available)
library(readr)
library(dplyr)

# use_r("topGO_All")
# use_r("load_topGO_db")
# use_r("RNAseqEasy")
# use_r("load_data")
# use_r("PCA_analysis")
#
# use_mit_license()
#
# use_testthat()
#
# use_test(name = "topGO_All")
# test_file("../tests/testthat/test-topGO_All.R")
# use_test(name = "load_data")
# use_test(name = "PCA_analysis")
#
# library(tidyverse)
#
# Mpo_GO_GOSLIM <- read_delim("G:/Mi unidad/Antiguo Drive/CNB/Marchantia/v6.1r1/Blazquez GO_db_no1.csv",
#                             delim = "\t", escape_double = FALSE,
#                             trim_ws = TRUE, col_names = TRUE)
#
# use_data(Mpo_GO_GOSLIM, overwrite = TRUE)

# DEG_example_list <- read.delim("G:/Mi unidad/Antiguo Drive/CNB/Ayuda Wendy/RNAseq/DESeq2/1. gh3aMock_vs_Tak1Mock/gh3aMock_vs_Tak1Mock_Sig.txt",
#                                header = TRUE)
# use_data(DEG_example_list, overwrite = TRUE)
#
# use_rmarkdown_template("topGO_All")
#
# use_vignette(name = "topGO_All", title = "Simplify topGO Enrichment and plotting")
#
#
# ## List of other packages to import when using this package
# use_package("dplyr", type = "Imports")
# use_package("tidyr", type = "Imports")
# use_package("ggplot2", type = "Imports")
# use_package("readr", type = "Imports")
# use_package("tximport", type = "Imports")
# use_package("DESeq2", type = "Imports")
# use_package("pheatmap", type = "Imports")
# use_package("topGO", type = "Imports")
# use_package("WGCNA", type = "Imports")
# use_package("GO.db", type = "Imports")
# use_package("rrvgo", type = "Imports")
# use_package("ggsci", type = "Imports")
# use_package("Hmisc", type = "Imports")
# use_package("openxlsx", type = "Imports")
# use_package("phylotools", type = "Imports")
#
#
#
# library(tximport)
# library(DESeq2)
# library(tidyverse)
# library(phylotools)
# library(pheatmap)
# library(topGO)
# library(WGCNA)
# library(Hmisc)
# library(openxlsx)
# library(wordcloud2)
# library(GO.db)
# library(htmlwidgets)
# library(rrvgo)
# library(pals)
# library(reshape2)
# library(ggsci)
# library(factoextra)
#
#

# sampleDir <- "/Users/salva/Google Drive/My Drive/Antiguo Drive/CNB/Ayuda Wendy/RNAseq/02_Salmon/"
# sample_table <- read.delim("/Users/salva/Google Drive/My Drive/Antiguo Drive/CNB/Ayuda Wendy/RNAseq/Sample_Data_Wendy.txt")
# Marchantia7_Transcripts <- get.fasta.name("/Users/salva/Google Drive/My Drive/Antiguo Drive/CNB/Marchantia/v7.1/MpTak_v7.1.mRNA.fa", clean_name = FALSE)
#
# use_data(Marchantia7_Transcripts, overwrite = TRUE)
#
# Marchantia7_tx2gene <- data.frame(Name = Marchantia7_Transcripts) %>%
#   separate(Name, sep = " ", into = c("TXNAME", "CDS")) %>%
#   mutate(GENEID = str_sub(TXNAME, 1, -3)) %>%
#   dplyr::select(-CDS)
#
# use_data(Marchantia7_tx2gene, overwrite = TRUE)
#
# Example_Wendy <- load_tximport_data(samplesDir = sampleDir, sample_table = sample_table, tx2gene = Marchantia7_tx2gene)
