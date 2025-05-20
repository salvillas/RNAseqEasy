library(usethis)
library(roxygen2)
library(devtools)
library(testthat)
library(available)

available("ggplot2")
usethis::use_git()

usethis::use_github()

use_r("topGO_All")
use_r("RNAseqEasy")

use_mit_license()

library(tidyverse)

Mpo_GO_GOSLIM <- read_delim("G:/Mi unidad/Antiguo Drive/CNB/Marchantia/v6.1r1/Blazquez GO_db_no1.csv",
                            delim = "\t", escape_double = FALSE,
                            trim_ws = TRUE, col_names = TRUE)

use_data(Mpo_GO_GOSLIM)

use_rmarkdown_template("topGO_All")


## List of other packages to import when using this package
use_package("dplyr", type = "Imports")
use_package("ggplot2", type = "Imports")
use_package("readr", type = "Imports")
use_package("tximport", type = "Imports")
use_package("DESeq2", type = "Imports")
use_package("pheatmap", type = "Imports")
use_package("topGO", type = "Imports")
use_package("WGCNA", type = "Imports")
use_package("GO.db", type = "Imports")
use_package("rrvgo", type = "Imports")
use_package("ggsci", type = "Imports")
use_package("Hmisc", type = "Imports")
use_package("openxlsx", type = "Imports")



library(tximport)
library(DESeq2)
library(tidyverse)
library(phylotools)
library(pheatmap)
library(topGO)
library(WGCNA)
library(Hmisc)
library(openxlsx)
library(wordcloud2)
library(GO.db)
library(htmlwidgets)
library(rrvgo)
library(pals)
library(reshape2)
library(ggsci)
library(factoextra)


