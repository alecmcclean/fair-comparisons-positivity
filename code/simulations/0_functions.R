###############################################################################
### Author: Alec McClean
### Purpose: Functions for smooth approximations
###############################################################################

# Propensity score function with varying margin condition, dictated by alpha
propensity_score <- function(x, alpha) {
  
  if (alpha == "infinity") {
    pi <- case_when(
      x < 0.1 ~ 0,
      x > 0.9 ~ 1,
      T ~ x
    )
  } else {
    constant <- 0.1
    beta <- 10^(1/alpha) * constant
    pi <- case_when(
      x < 0.1 ~ 0, # Set to zero below 0.1
      between(x, 0.1, 0.2) ~ beta * (x - 0.1)^(1/alpha), # Impose margin between 0.1 and 0.2
      between(x, 0.2, 0.8) ~ constant + (x - 0.2) * (4 / 3), # Linearly interpolate from 0.2 to 0.8
      between(x, 0.8, 0.9) ~ 1 - beta * (0.9 - x)^(1/alpha), # Impose margin condition between 0.8 and 0.9
      T ~ 1 # Set to one above 0.9
    )
  }
  return(pi)
}

# CATE
CATE <- function(x) { x + x^2 }

# Smooth approximation of indicator
smooth_exp <- function(x, k) 1 - exp(-1 * k * x)
smooth_exp_diff <- function(x, k) k * exp(-1 * k * x)
smooth_approx <- function(x, k) smooth_exp(x, k) * smooth_exp(1 - x, k)


# Construct STATE estimator's IF values
state_if_values <- function(smooth_func, smooth_func_deriv, A, Y, pihat, muhat0, muhat1) {
  
  S <- smooth_func(pihat) * smooth_func(1 - pihat)
  varphi_S <- smooth_func_deriv(1 - pihat) * (pihat - A) * smooth_func(pihat) +
    smooth_func_deriv(pihat) * (A - pihat) * smooth_func(1 - pihat)
  
  varphi_mu1 <- ifelse(S == 0, 0, (S * A / pihat) * (Y - muhat1))
  varphi_mu0 <- ifelse(S == 0, 0, (S * (1 - A) / (1 - pihat)) * (Y - muhat0))
  varphi_ate_S <- varphi_mu1 - varphi_mu0
  
  if_values <- (muhat1 - muhat0) * (S + varphi_S) + varphi_ate_S
  return(if_values)
  
}

construct_data <- function(N, ALPHA) {
  
  # Construct data
  X <- runif(n = N)
  pi <- propensity_score(X, alpha = ALPHA)
  A <- rbinom(n = N, size = 1, prob = pi)
  Y0 <- rnorm(n = N) # Set E(Y^0 | X) = 0 
  Y1 <- CATE(X) + rnorm(n = N) # E(Y^1 - Y^0 | X) = X + X^2 (i.e., heterogeneity)
  Y <- A * Y1 + (1 - A) * Y0
  
  # Construct data frame
  dat <- data.frame(
    X = X, pi = pi, A = A, Y = Y, Y0 = Y0, Y1 = Y1
  )
  
  return(dat)
}

run_simulation <- function(K, N, ALPHA, RATE) {
  
  if (ALPHA == "infinity") {
    ALPHA <- "infinity"
  } else {
    ALPHA <- as.numeric(as.character(ALPHA))
  }
  
  # Construct data
  dat <- construct_data(N, ALPHA)
  
  # Compute STATE
  if (ALPHA == "infinity") {
    STATE <- integrate(function(X) CATE(X), lower = 0.1, upper = 0.9)$value
  } else {
    STATE <- integrate(
      function(X) CATE(X) * smooth_approx(x = propensity_score(x = X, alpha = ALPHA),
                                          k = K),
      lower = 0, upper = 1
    )$value
  }
  
  if (RATE == 0) {
    
    # When rate == 0, assume no error.
    pihat <- dat$pi
    muhat0 <- 0
    muhat1 <- CATE(dat$X)
    
  } else {
    
    # Construct propensity score estimator
    direction_error <- sample(c(-1, 1), N, replace = TRUE)
    pihat <- dat$pi + direction_error * 10^(-1) * rnorm(n = N, mean = N^(-1 * RATE), sd = N^(-1 * RATE))
    pihat <- pmin(1, pmax(0, pihat))
    
    if (ALPHA == "infinity") {
      pihat <- ifelse(pihat < 0.1, 0,
                      ifelse(pihat > 0.9, 1, pihat))
    }
    
    # Construct regression function estimators
    direction_error <- sample(c(-1, 1), N, replace = TRUE)
    muhat0 <- 0 + direction_error * 10^(-1) * rnorm(n = N, mean = N^(-1 * RATE), sd = N^(-1 * RATE))
    
    direction_error <- sample(c(-1, 1), N, replace = TRUE)
    muhat1 <- CATE(dat$X) + direction_error * 10^(-1) * rnorm(n = N, mean = N^(-1 * RATE), sd = N^(-1 * RATE))
    
  }
  
  # Compute influence function values
  state_ifvalues <- state_if_values(
    smooth_func = function(x) smooth_exp(k = K, x),
    smooth_func_deriv = function(x) smooth_exp_diff(k = K, x),
    A = dat$A, Y = dat$Y, pihat = pihat, muhat0 = muhat0, muhat1 = muhat1
  )
  
  # Compile results: point estimate and confidence intervals
  results <- data.frame(
    k = K, n = N, state = STATE, alpha = as.character(ALPHA), rate = RATE,
    state_est = mean(state_ifvalues), 
    smooth_tate_est = mean(state_ifvalues),
    state_ci_ll = mean(state_ifvalues) - qnorm(0.975) * sd(state_ifvalues) / sqrt(N),
    state_ci_ul = mean(state_ifvalues) + qnorm(0.975) * sd(state_ifvalues) / sqrt(N)
  )
  
  return(results)
}
