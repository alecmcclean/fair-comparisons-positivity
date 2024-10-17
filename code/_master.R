###############################################################################
### Author: Alec McClean
### Purpose: Master script for fair rankings of providers
###############################################################################

if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(magrittr,
       tidyverse,
       ggthemes,
       latex2exp)

options(stringsAsFactors = F)

set.seed(20240903)

source("1_analysis.R")
source("2_figures.R")