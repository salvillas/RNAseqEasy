library(usethis)
library(roxygen2)
library(devtools)
library(testthat)
usethis::use_git()

usethis::use_github()

use_r("topGO_All")
use_r("RNAseqEasy")

library(tidyverse)

Mpo_GO_GOSLIM <- read_delim("G:/Mi unidad/Antiguo Drive/CNB/Marchantia/v6.1r1/Blazquez GO_db_no1.csv",
                            delim = "\t", escape_double = FALSE,
                            trim_ws = TRUE, col_names = TRUE)

use_data(Mpo_GO_GOSLIM)

use_rmarkdown_template("topGO_All")
