###############################################################################
### Author: Alec McClean
### Purpose: Make figures
###############################################################################

ny_data_top10 <- readRDS("../data/ny_data_top10.RDS")
load("../output/utils.RData")

# Positivity violations
plot_data <- nuis %>% 
  as.data.frame %>% 
  select(contains("pi")) %>% 
  pivot_longer(cols = everything(),
               names_to = "provider",
               values_to = "pihat") %>%
  mutate(provider = gsub("pi_a", "", provider)) 

# Number of observations with near or absolute positivity violations
mean(plot_data$pihat < 0.01); sum(plot_data$pihat == 0)

# Result figure
files <- paste0("../output/estimates_", c("shift_d0.90", "shift_d0.50", 
                                          "tilt_d0.90",  "tilt_d0.50",
                                          "tsm"), "_k10^2.RDS")
read_add_filename <- function(file) {
  data <- readRDS(file)
  
  data$res %>%
    mutate(filename = file)
}

combined_data <- lapply(files, read_add_filename) %>% bind_rows()
change_name <- data.frame(filename = files, 
                          spec = c("Shift: delta = 0.9",
                                  "Shift: delta = 0.5",
                                  "Tilt: delta = 0.9",
                                  "Tilt: delta = 0.5",
                                  "Smooth trimmed TSMs"))

combined_data %<>% full_join(change_name) %>% select(-filename) 
combined_data %<>% mutate(spec = factor(spec, levels = c("Shift: delta = 0.9",
                                                         "Shift: delta = 0.5",
                                                         "Tilt: delta = 0.9",
                                                         "Tilt: delta = 0.5",
                                                         "Smooth trimmed TSMs")))

est_noint <- readRDS("../output/estimates_noint.RDS")
est_noint <- est_noint$res$est[1]

### Anonymize the ID
anonymized_id <- ny_data_top10 %>% select(PROVUSRD) %>% unique() %>% 
  mutate(PROVUSRD = as.character(PROVUSRD))
anonymized_id$anon_id <- as.roman(1:10)
combined_data %<>% left_join(anonymized_id) %>% select(-PROVUSRD)

combined_data %<>% filter(spec == "Smooth trimmed TSMs")

results_plot <-
  ggplot(data = combined_data, aes(x = reorder(as.character(anon_id), anon_id), y = est)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci.ll, ymax = ci.ul)) +
  geom_hline(yintercept = est_noint, linetype = "dashed", color = "red") +
  theme_minimal() +
  theme(text = element_text(family = "serif")) +
  labs(x = "Anonymized provider ID", 
       y = "Point estimate with pointwise 95% confidence interval")

ggsave(plot = results_plot, filename = "../figures/results.png",
       width = 5, height = 5, bg = "white")

################################
### Differences for trimmed TSMs

trimmed_results <- readRDS("../output/estimates_tsm_k10^2.RDS")
cov_mat <- trimmed_results$covmat
res <- trimmed_results$res

# Focus on the difference between providers with max and minimum effect
max_index <- which(res$est == max(res$est))
full_index <- seq(1, length(res$PROVUSRD), 1) 
full_index <- full_index[full_index != max_index]

for (other_index in full_index) {
  diff <- max(res$est) - res$est[other_index]
  se <- sqrt(cov_mat[max_index, max_index] + cov_mat[other_index, other_index] -
               2 * cov_mat[max_index, other_index])
  
  ci_ll <- diff - qnorm(0.975) * se
  ci_ul <- diff + qnorm(0.975) * se
  
  cat("\nEstimate of provider ", 
      as.character(anonymized_id$anon_id[anonymized_id$PROVUSRD == res$PROVUSRD[max_index]]), " readmission rate minus provider ",
      as.character(anonymized_id$anon_id[anonymized_id$PROVUSRD == res$PROVUSRD[other_index]]), " readmission rate: ", diff,
      "\n95% CI: [", ci_ll, ", ", ci_ul, " ]")
}

rm(list = ls(all = T))
gc()


