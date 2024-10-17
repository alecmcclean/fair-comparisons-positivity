###############################################################################
### Author: Alec McClean
### Purpose: Estimate fair ranking of effects.
###############################################################################

source("0_functions.R")

############################
### Load and clean data
############################

load("../data/USRDS_2020_inpatient.RData")
ny_data <- filter(data, STATE == "36")

# Focus on just 50 largest providers
top_hospitals <- ny_data %>%
  count(PROVUSRD) %>%
  arrange(desc(n)) %>%
  slice_head(n = 10) %>%
  pull(PROVUSRD)

ny_data_top10 <- ny_data %>%
  filter(PROVUSRD %in% top_hospitals)

saveRDS(ny_data_top10, file = "../data/ny_data_top10.RDS")
rm(ny_data, data)


############################
### Run estimators
############################

# Specify constant in smooth approximation
K <- 10^2

# Estimate shifted effects: multiplicative with $\delta = 0.9$
est_shift <- tsfx(
  y = ny_data_top10$readmit30,
  a = ny_data_top10$PROVUSRD,
  x = ny_data_top10 %>% select(-readmit30, -PROVUSRD, -STATE, -CLM_THRU),
  splits = 2,
  sl.lib = c("SL.mean", "SL.ranger", "SL.gam", "SL.glm", "SL.glm.interaction"),
  smooth_func = function(x) smooth_exp(x, K),
  smooth_func_diff = function(x) K * exp(-1 * K * x),
  int_func = function(x) shift(x, delta = 0.9),
  int_func_diff = function(x) shift_diff(x, delta = 0.9))

# Save everything
saveRDS(est_shift, "../output/estimates_shift_d0.90_k10^2.RDS")
nuis <- est_shift[["nuis"]]; utils <- est_shift[["utils"]]
save(nuis, utils, file = "../output/utils_n50.RData")

# Re-use outputs from previous estimation but with different shift functions
# Multiplicative shift with $\delta = 0.5$
est_shift2 <- reanalyze_output(nuis, utils,
                               smooth_func = function(x) smooth_exp(x, K),
                               smooth_func_diff = function(x) K * exp(-1 * K * x),
                               int_func = function(x) shift(x, delta = 0.5),
                               int_func_diff = function(x) shift_diff(x, delta = 0.5))

saveRDS(est_shift2, "../output/estimates_shift_d0.50_k10^2.RDS")

# Exponential tilt with $\delta = 0.9$
est_tilt <- reanalyze_output(nuis, utils,
                             smooth_func = function(x) smooth_exp(x, K),
                             smooth_func_diff = function(x) K * exp(-1 * K * x),
                             int_func = function(x) tilt(x, delta = 0.9),
                             int_func_diff = function(x) tilt_diff(x, delta = 0.9))

saveRDS(est_tilt, "../output/estimates_tilt_d0.90_k10^2.RDS")

# Exponential tilt with $\delta = 0.5$
est_tilt2 <- reanalyze_output(nuis, utils,
                             smooth_func = function(x) smooth_exp(x, K),
                             smooth_func_diff = function(x) K * exp(-1 * K * x),
                             int_func = function(x) tilt(x, delta = 0.5),
                             int_func_diff = function(x) tilt_diff(x, delta = 0.5))

saveRDS(est_tilt2, "../output/estimates_tilt_d0.50_k10^2.RDS")

# Treatment-specific means
est_tsm <- reanalyze_output(nuis, utils,
                            smooth_func = function(x) smooth_exp(x, K),
                            smooth_func_diff = function(x) K * exp(-1 * K * x),
                            int_func = function(x) tsm(x),
                            int_func_diff = function(x) tsm_diff(x))

saveRDS(est_tsm, "../output/estimates_tsm_k10^2.RDS")

# No intervention
est_noint <- reanalyze_output(nuis, utils,
                            smooth_func = function(x) smooth_exp(x, K),
                            smooth_func_diff = function(x) K * exp(-1 * K * x),
                            int_func = function(x) shift(x, delta = 1),
                            int_func_diff = function(x) shift_diff(x, delta = 1))

saveRDS(est_noint, "../output/estimates_noint.RDS")

rm(list = ls(all = T))
gc()
