###############################################################################
### Author: Alec McClean
### Purpose: Simulations for estimating STATE
###############################################################################

set.seed(20250713)
source("0_functions.R")

#######################
### Example data

# Construct data
N <- 1000 # Sample size
ALPHA <- 0.5
K <- 20
X <- runif(n = N)
pi <- propensity_score(X, alpha = ALPHA)
A <- rbinom(n = N, size = 1, prob = pi)
Y0 <- rnorm(n = N) # Set E(Y^0 | X) = X 
Y1 <- CATE(X) + rnorm(n = N) # E(Y^1 - Y^0 | X) = X^2 (i.e., heterogeneity)
Y <- A * Y1 + (1-A) * Y0

dat <- data.frame(
  X = X, pi = pi, A = A, Y = Y, Y0 = Y0, Y1 = Y1
)

# Plot propensity score density and CATE
p1 <- ggplot(data = dat, aes(x = X, y = pi)) + 
  geom_line() +
  labs(x = "Covariate value", 
       y = "Propensity score") +
  theme_clean() +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")

ggsave(plot = p1, filename = "../../figures/simulations/positivity_violations.png",
       width = 6, height = 4, bg = "white")

p2 <- ggplot(data = dat, aes(x = X, y = X^2)) + geom_line() + 
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "red") + 
  geom_vline(xintercept = 0.9, linetype = "dashed", color = "red") +
  theme_clean() +
  labs(x = "Covariate value",
       y = "Conditional Average Treatment Effect")

ggsave(plot = p2, filename = "../../figures/simulations/cate.png",
       width = 6, height = 4, bg = "white")

dat %<>% mutate(S = smooth_exp(pi, k = K) * smooth_exp(1 - pi, k = K),
                ind = as.integer(0 < pi & pi < 1))

p3 <- dat %>% gather(var, value, ind, S) %>%
  mutate(var = ifelse(var == "ind", "Indicator", "Smooth indicator")) %>%
  ggplot(aes(x = X, y = value, color = var, group = var)) +
  geom_line() +
  theme_clean() +
  labs(color = "Trimming set", 
       x = "Covariate value",
       y = "Trimming indicator or smooth trimming weight")

ggsave(plot = p3, filename = "../../figures/simulations/trimming_indicators.png",
       width = 6, height = 4, bg = "white")

p <- (p1 + p2) / p3
ggsave(plot = p, filename = "../../figures/simulations/all_plots.png")


# Compute TATE and STATE
STATE <- integrate(
  function(X) CATE(X) * smooth_approx(x = propensity_score(x = X, alpha = ALPHA),
                                      k = K),
  lower = 0, upper = 1
)$value

cat("Smooth Trimmed ATE ( k = ", K, ", alpha = ", ALPHA, "): ", STATE)

# Clean
rm(list = ls(all = T))
gc()


###############################
### Examine STATE coverage

source("0_functions.R")
NUM_ITERS <- 200
specs <- expand.grid(iter = 1:NUM_ITERS,
                     N = 10^seq(2, 5, 1), 
                     K = 20,
                     ALPHA = 0.5,
                     RATE = c(0.2, 0.3, 0.4, 0.5))

results_dat <- data.frame()
for (row in 1:nrow(specs)) {
  cat("\nRow:", format(row, big.mark = ","), "out of", format(nrow(specs), big.mark = ","))
  results_dat %<>% bind_rows(
    bind_cols(
      run_simulation(K = specs$K[row],
                     N = specs$N[row],
                     ALPHA = specs$ALPHA[row],
                     RATE = specs$RATE[row]),
      data.frame(iter = specs$iter[row])
    )
  )
}

saveRDS(results_dat, "../../output/simulations/results_state_20250622.RDS")

# Calculate coverage
coverage <- results_dat %>% 
  group_by(k, n, alpha, rate) %>%
  summarize(
    coverage = mean(between(state, state_ci_ll, state_ci_ul)),
    coverage_ci_lb = coverage - qnorm(0.975) * sqrt((coverage * (1 - coverage)) / n()),
    coverage_ci_ub = coverage + qnorm(0.975) * sqrt((coverage * (1 - coverage)) / n())
  ) %>% 
  mutate(rate = paste0("Nuisance convergence rate: n^{-", rate, "}"))

p <- ggplot(data = coverage, aes(x = factor(n), y = coverage)) +
  geom_point() +
  geom_errorbar(aes(ymin = coverage_ci_lb, ymax = coverage_ci_ub)) +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed") +
  facet_wrap(~ rate) +
  theme_clean() +
  labs(x = "Sample size", y = "Coverage")


ggsave(p, filename = "../../figures/simulations/coverage.png",
       width = 8, height = 6, bg = "white")
