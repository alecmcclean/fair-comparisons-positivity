###############################################################################
### Author: Alec McClean
### Purpose: Functions for fair provider profiling
###############################################################################

# Smooth approximation of indicator
smooth_exp <- function(x, k) 1 - exp(-1 * k * x)
smooth_exp_diff <- function(x, k) k * exp(-1 * k * x)

# Tilt, shift, and TSM, and their derivatives
tilt <- function(x, delta) (delta * x) / (delta * x + 1 - x)
tilt_diff <- function(x, delta) delta * (1 - x) / (delta * x + 1 - x)^2

shift <- function(x, delta) delta * x
shift_diff <- function(x, delta) delta

tsm <- function(x) 0
tsm_diff <- function(x) 0

tsfx <- function(y, a, x, splits = 2,
                 sl.lib = c("SL.earth", "SL.gam", "SL.glm", "SL.glm.interaction",
                            "SL.mean", "SL.ranger", "SL.rpart"),
                 smooth_func,  # Smoothing function for trimming
                 smooth_func_diff, # Derivative of smoothing function
                 int_func, # Interventional function
                 int_func_diff) {   # Interventional function derivative
  
  # Load required libraries
  require("SuperLearner")
  require("earth")
  require("gam")
  require("ranger")
  require("rpart")

  ################################
  ### Setup
  ################################
  
  # Set up some basic variables
  n <- dim(x)[1]  # Number of observations
  avals <- names(table(a))  # Unique treatment levels
  n.avals <- length(avals)  # Number of unique treatments
  pb <- txtProgressBar(min = 0, max = 2 * splits * n.avals, style = 3)  # Progress bar
  
  # Create cross-validation splits
  folds <- sample(rep(1:splits, ceiling(n / splits))[1:n])
  
  # Initialize matrices to store predicted outcomes (muhat) and propensity scores (pihat)
  muhat <- as.data.frame(matrix(NA, nrow = n, ncol = n.avals))
  colnames(muhat) <- paste("a", avals, sep = "")
  pihat <- muhat
  
  ################################
  ### Estimate nuisance functions
  ################################
  
  pbcount <- 0  # Initialize progress bar counter
  for (i in 1:n.avals) {  # Loop over each treatment level
    if (i == 1) { Sys.sleep(0.1); setTxtProgressBar(pb, pbcount); pbcount <- pbcount + 1 }
    for (vfold in 1:splits) {  # Loop over each cross-validation fold
      
      # Create training and test sets
      train <- folds != vfold; test <- folds == vfold
      if (splits == 1) { train <- test }  # If no splits, use all data for training
      
      # Estimate propensity score, but not for the last treatment level
      # For last treatment level, see below
      pifit <- SuperLearner(as.numeric(a == avals[i])[train], as.data.frame(x[train,]),
                            newX = as.data.frame(x[test,]), SL.library = sl.lib, family = binomial)
      pihat[test, i] <- pifit$SL.predict
    
      
      Sys.sleep(0.1)
      setTxtProgressBar(pb, pbcount); pbcount <- pbcount + 1
      
      # Estimate regression function for the current treatment level
      mufit <- SuperLearner(y[a == avals[i] & train],
                            as.data.frame(x[a == avals[i] & train,]),
                            newX = as.data.frame(x[test,]), SL.library = sl.lib)
      muhat[test, i] <- mufit$SL.predict
      
      Sys.sleep(0.1)
      setTxtProgressBar(pb, pbcount); pbcount <- pbcount + 1
    }
  }
  
  # Normalize propensity score values to sum to 1.
  pihat <- sweep(pihat, 1, rowSums(pihat), FUN = "/")

  #############################################
  ### Calculate influence function values
  #############################################

  ifvalues <- calculate_if_values(pihat, a, avals, n.avals, muhat, y, n,
                                  smooth_func, smooth_func_diff,
                                  int_func, int_func_diff)
    
  # Compute estimates, standard errors, confidence intervals, and p-values
  tryCatch({
    est <- apply(ifvalues[["uncentered_ifmat"]], 2, mean)
  }, error = function(e) {
      # Enter debugging mode if an error occurs
      message("An error occurred: ", e$message)
      browser()  # Enter browser mode for debugging
    })
  
  se <- apply(ifvalues[["uncentered_ifmat"]], 2, sd) / sqrt(n)
  ci.ll <- est - qnorm(0.975) * se; ci.ul <- est + qnorm(0.975) * se
  res <- data.frame(PROVUSRD = avals, est, se, ci.ll, ci.ul)

  # Covariance matrix. Can be used for joint confidence intervals
  cov_mat <- cov(ifvalues[["centered_ifmat"]]) / n
  
  Sys.sleep(0.1)
  setTxtProgressBar(pb, pbcount)
  close(pb)
  
  # Combine nuisance parameter estimates and return results invisibly
  nuis <- as.data.frame(cbind(pihat, muhat))
  colnames(nuis) <- paste(rep(c("pi", "mu"), rep(n.avals, 2)), colnames(nuis), sep = "_")
  
  # Other utilities
  utils <- list(avals = avals, n.avals = n.avals, a = a, y = y, n = n)
  
  print(res)
  return(invisible(list(res = res, covmat = cov_mat, nuis = nuis, utils = utils)))
}

calculate_if_values <- function(pihat, a, avals, n.avals, 
                                muhat, y, n,
                                smooth_func, smooth_func_diff,
                                int_func, int_func_diff) {
  
  # Calculate smoothed propensity scores
  smooth_pihat <- smooth_func(pihat)
  smooth_pihat_diff <- smooth_func_diff(pihat)
  
  # Calculate S(X \in C_X)
  S <- apply(smooth_pihat, 1, prod)
  
  # Calculate EIF of S(X \in C_X)
  EIF_S <- rep(NA, n)
  for (subject in 1:n) {
    summat <- matrix(rep(pihat[subject, ], n.avals), nrow = n.avals, byrow = TRUE)
    for (txnum in 1:n.avals) {
      summat[txnum, txnum] <- smooth_func_diff(pihat[subject, txnum]) * 
        (as.integer(a[subject] == avals[txnum]) - pihat[subject, txnum])
    }
    EIF_S[subject] <- 
      sum(apply(matrix(as.numeric(summat), nrow = nrow(summat), byrow = TRUE), 1, prod))
  }
  
  # Calculate Q_a(A=b \mid X) and its EIF
  q_a <- array(1:(n * n.avals * n.avals), dim = c(n, n.avals, n.avals))
  eif_q_a <- Q_a <- eif_Q_a <- q_a
  
  for (subject in 1:n) { # For each subject
    for (txnum_a in 1:n.avals) { # Consider intervention targeting a
      for (txnum_b in 1:n.avals) { # Calculate interventional propensity score at b
        tx_a <- avals[txnum_a]
        tx_b <- avals[txnum_b]
        q_a[subject, txnum_a, txnum_b] <- 
          ifelse(tx_b != tx_a,
                 int_func(pihat[subject, txnum_b]),
                 1 - sum(int_func(pihat[subject, -txnum_a])))
        
        eif_q_a[subject, txnum_a, txnum_b] <- 
          ifelse(tx_b != tx_a,
                 int_func_diff(pihat[subject, txnum_b]) *
                   (as.integer(a[subject] == tx_b) - pihat[subject, txnum_b]),
                 -1 * sum(int_func_diff(pihat[subject, -txnum_a]) * 
                            (as.integer(a[subject] == avals[-txnum_a]) - pihat[subject, -txnum_a])))
        
        # Calculate interventional propensity score
        Q_a[subject, txnum_a, txnum_b] <- 
          S[subject] * q_a[subject, txnum_a, txnum_b] + (1 - S[subject]) * pihat[subject, txnum_b]
        
        # Calculate EIF of interventional propensity score
        eif_Q_a[subject, txnum_a, txnum_b] <-
          (1 - S[subject]) * (as.integer(a[subject] == tx_b) - pihat[subject, txnum_b]) -
          EIF_S[subject] * pihat[subject, txnum_b] +
          EIF_S[subject] * q_a[subject, txnum_a, txnum_b] + 
          S[subject] * eif_q_a[subject, txnum_a, txnum_b]
      }
    }
  } 
  
  # Calculate influence function values for each intervention
  centered_ifmat <- matrix(rep(a, n.avals), nrow = n, byrow = FALSE)
  uncentered_ifmat <- centered_ifmat
  for(subject in 1:n) {
    for(txnum_a in 1:n.avals) {
      
      # Q / pi. If we have positive number divided by zero, return NA. Otherwise,
      # set 0 / 0 = 0 or return Q / pi.
      resid_weight <- 
        unlist(ifelse(Q_a[subject, txnum_a, ] > 0 & pihat[subject, ] == 0,
                      NA,
                      ifelse(Q_a[subject, txnum_a, ] == 0, 0,
                             Q_a[subject, txnum_a, ] / pihat[subject, ])))
      
      # Stop if we created any NAs. The effects are not robust to positivity violations!
      stopifnot(!is.na(resid_weight))
      
      tryCatch({
        centered_ifmat[subject, txnum_a] <- sum(
          resid_weight * as.integer(a[subject] == avals) * (y[subject] - muhat[subject, ]) + 
            eif_Q_a[subject, txnum_a, ] * muhat[subject, ]
        )}, error = function(e) {
          # Enter debugging mode if an error occurs
          message("An error occurred: ", e$message)
          browser()  # Enter browser mode for debugging
        })
      
      uncentered_ifmat[subject, txnum_a] <- centered_ifmat[subject, txnum_a] +
        sum(Q_a[subject, txnum_a, ] * muhat[subject, ]) 
      
    }
  }
  
  return(list(centered_ifmat = centered_ifmat, uncentered_ifmat = uncentered_ifmat))
}

reanalyze_output <- function(nuis, utils,
                             smooth_func, smooth_func_diff,
                             int_func, int_func_diff) {
  
  pihat <- nuis %>% select(contains("pi")) %>% as.matrix
  muhat <- nuis %>% select(contains("mu")) %>% as.matrix
  
  avals <- utils$avals
  n.avals <- utils$n.avals
  a <- utils$a
  y <- utils$y
  n <- utils$n
  
  ifvalues <- calculate_if_values(pihat, a, avals, n.avals, muhat, y, n,
                                  smooth_func, smooth_func_diff,
                                  int_func, int_func_diff)
  
  # Compute estimates, standard errors, confidence intervals, and p-values
  tryCatch({
    est <- apply(ifvalues[["uncentered_ifmat"]], 2, mean)
  }, error = function(e) {
    # Enter debugging mode if an error occurs
    message("An error occurred: ", e$message)
    browser()  # Enter browser mode for debugging
  })
  
  se <- apply(ifvalues[["uncentered_ifmat"]], 2, sd) / sqrt(n)
  ci.ll <- est - qnorm(0.975) * se; ci.ul <- est + qnorm(0.975) * se
  res <- data.frame(PROVUSRD = avals, est, se, ci.ll, ci.ul)
  
  # Covariance matrix. Can be used for joint confidence intervals
  cov_mat <- cov(ifvalues[["centered_ifmat"]]) / n

  # Combine nuisance parameter estimates and return results invisibly
  print(res)
  return(invisible(list(res = res, covmat = cov_mat)))
}
