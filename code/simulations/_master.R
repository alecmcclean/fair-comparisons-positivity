###############################################################################
### Author: Alec McClean
### Purpose: Master script for smoothed trimmed ATE
###############################################################################

if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(magrittr,
       tidyverse,
       ggthemes,
       latex2exp,
       patchwork)

options(stringsAsFactors = F)

source("1_simulations.R")