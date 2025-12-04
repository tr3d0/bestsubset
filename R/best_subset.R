#' Best Subset Selection via MIQP (Fixed k)
#'
#' @description
#' Solves the best subset selection problem for a fixed subset size \code{k} using
#' Mixed-Integer Quadratic Programming (MIQP).
#'
#' @details
#' This function formulates the variable selection problem as a quadratic optimization
#' problem with linear constraints. It leverages the 'Gurobi' optimizer to find the
#' globally optimal subset of predictors that minimizes the Residual Sum of Squares (RSS).
#' 
#' The method handles high-dimensional data by using a Big-M formulation to constrain
#' the regression coefficients.
#'
#' @param X Numeric matrix. The design matrix (predictors).
#' @param y Numeric vector. The response vector.
#' @param k Integer. The exact number of variables to select.
#' @param M Numeric. Big-M parameter for constraints. Default is 15.
#' @param time_limit Numeric. Time limit for the solver in seconds. Default is 300.
#' @param mip_gap Numeric. The relative MIP optimality gap. Default is 0.01 (1%).
#' @param verbose Logical. If TRUE, shows Gurobi console output. Default is TRUE.
#' @param threads Integer. Number of threads to use (0 = all available). Default is 0.
#' @param mip_focus Integer. Gurobi MIPFocus parameter (1=Feasibility, 2=Optimality, 3=Bound).
#' @param presolve Integer. Gurobi Presolve level (-1=Auto, 0=Off, 1=Conservative, 2=Aggressive).
#' @param ... Additional parameters passed directly to the Gurobi solver.
#'
#' @return A list containing:
#' \item{coefficients}{Estimated coefficients (including intercept).}
#' \item{beta_vector}{Estimated coefficients (excluding intercept).}
#' \item{selected_indices}{Indices of the selected variables.}
#' \item{rss}{Residual Sum of Squares on original scale.}
#' \item{gap}{MIP gap achieved.}
#' \item{status}{Gurobi optimization status.}
#'
#' @importFrom Matrix crossprod Diagonal Matrix as
#' @importFrom stats start
#' @import gurobi
#' @export
best_subset <- function(X, y, k, M = 15, 
                                  time_limit = 300, 
                                  mip_gap = 0.01,
                                  verbose = TRUE,
                                  threads = 0,
                                  mip_focus = 3,
                                  presolve = 2,
                                  ...) {

  # 1. Data Preparation (Standardization)
  X_scaled <- scale(X, center = TRUE, scale = TRUE)
  y_scaled <- scale(y, center = TRUE, scale = FALSE)

  scale_X <- attr(X_scaled, 'scaled:scale')
  center_X <- attr(X_scaled, 'scaled:center')
  mean_y <- mean(y)

  n <- nrow(X_scaled)
  p <- ncol(X_scaled)

  # 2. Quadratic Objective
  XtX <- crossprod(X_scaled)
  Xty <- crossprod(X_scaled, y_scaled)
  XtX <- as(XtX, "generalMatrix")

  # Q matrix: [2X'X, 0; 0, 0]
  Q <- Matrix(0, 2 * p, 2 * p, sparse = TRUE)
  Q[1:p, 1:p] <- 2 * XtX

  # Linear term: -2X'y for beta, 0 for z
  obj <- c(-2 * Xty, rep(0, p))

  # 3. Constraints (Big-M)
  # beta - Mz <= 0
  A1 <- cbind(Diagonal(p), Diagonal(p, x = -M))
  # -beta - Mz <= 0
  A2 <- cbind(Diagonal(p, x = -1), Diagonal(p, x = -M))
  # sum(z) = k
  A3 <- c(rep(0, p), rep(1, p))

  A <- rbind(A1, A2, A3)
  rhs <- c(rep(0, 2 * p), k)
  sense <- c(rep("<", 2 * p), "=")

  # 4. Gurobi Model
  model <- list()
  model$Q <- Q
  model$obj <- obj
  model$A <- A
  model$rhs <- rhs
  model$sense <- sense
  model$vtype <- c(rep("C", p), rep("B", p))

  # 5. Parameters Setup
  # Base parameters
  params <- list(
    OutputFlag = as.numeric(verbose),
    TimeLimit = time_limit,
    MIPGap = mip_gap,
    MIPFocus = mip_focus,
    Presolve = presolve,
    Threads = threads
  )
  
  # Append extra parameters passed via ... (if any)
  extra_params <- list(...)
  if (length(extra_params) > 0) {
    params <- c(params, extra_params)
  }

  # 6. Solve
  result <- gurobi(model, params)

  # 7. Process Results
  if (result$status %in% c("OPTIMAL", "TIME_LIMIT", "SUBOPTIMAL")) {
    beta_hat_scaled <- result$x[1:p]
    z_hat <- result$x[(p + 1):(2 * p)]

    selected_vars <- which(z_hat > 0.5)

    # Back-scaling
    beta_final <- rep(0, p)
    if(length(selected_vars) > 0) {
      beta_final[selected_vars] <- beta_hat_scaled[selected_vars] / scale_X[selected_vars]
    }

    intercept <- mean_y - sum(center_X * beta_final)

    # RSS on original scale
    y_pred <- intercept + X %*% beta_final
    rss_real <- sum((y - y_pred)^2)

    return(list(
      coefficients = c(Intercept = intercept, beta_final),
      beta_vector = beta_final,
      selected_indices = selected_vars,
      rss = rss_real,
      gap = result$mipgap,
      status = result$status
    ))
  } else {
    return(list(status = result$status, rss = Inf, selected_indices = integer(0)))
  }
}

#' Best Subset Selection with Automatic k Selection (EBIC)
#'
#' @description
#' Runs the best subset selection algorithm for a range of subset sizes (from 1 to \code{k_max})
#' and selects the globally optimal model minimizing the Extended Bayesian Information Criterion (EBIC).
#'
#' @details
#' This function acts as a wrapper around \code{\link{best_subset}}. It iterates through
#' all possible subset sizes up to \code{k_max}. For each size \code{k}, it fits the best subset
#' model using the MIQP solver.
#'
#' After fitting, it calculates the EBIC for each model using the formula from Chen & Chen (2008):
#' \deqn{EBIC = n \log(RSS/n) + k \log(n) + 2 \gamma \ln \binom{p}{k}}
#' where \eqn{\gamma = 1} if \eqn{p > n} (high-dimensional settings) and \eqn{\gamma = 0} otherwise.
#'
#' The function returns the model corresponding to the minimum EBIC value.
#'
#' @param X Numeric matrix. The design matrix (predictors).
#' @param y Numeric vector. The response vector.
#' @param k_max Integer. Maximum subset size to consider. Default is 10.
#' @param time_limit_per_k Numeric. Time limit per k in seconds. Default is 100.
#' @param mip_gap Numeric. The relative MIP optimality gap. Default is 0.01.
#' @param verbose Logical. If TRUE, shows Gurobi output for EACH step (can be very noisy). Default is FALSE.
#' @param ... Additional parameters passed to \code{\link{best_subset}} (e.g., \code{threads}, \code{presolve}).
#'
#' @return A list containing the best model's details (as returned by \code{\link{best_subset}}) augmented with:
#' \item{best_k_selected}{The integer k that minimized EBIC.}
#' \item{all_ebic}{A vector of EBIC scores for each k.}
#' \item{all_rss}{A vector of RSS values for each k.}
#'
#' @export
best_subset_auto <- function(X, y, k_max = 10, 
                             time_limit_per_k = 100, 
                             mip_gap = 0.01, 
                             verbose = FALSE, 
                             ...) {

  n <- nrow(X)
  p <- ncol(X)

  results_list <- list()
  ebic_scores <- numeric(k_max)
  rss_scores <- numeric(k_max)

  # EBIC gamma parameter (Chen \& Chen 2008)
  gamma_ebic <- ifelse(p > n, 1, 0)

  cat("Starting Best Subset Selection (k_max =", k_max, ")...\n")

  for (k in 1:k_max) {
    cat(sprintf("  Solving for k = %d ... ", k))

    # Pass parameters down to the optimizer
    # Note: explicit args overwrite defaults, ... passes the rest
    fit <- best_subset(X, y, k = k, 
                                 time_limit = time_limit_per_k, 
                                 mip_gap = mip_gap,
                                 verbose = verbose,
                                 ...)

    if (fit$status %in% c("OPTIMAL", "TIME_LIMIT", "SUBOPTIMAL")) {
      rss <- fit$rss
      if (rss <= 1e-10) rss <- 1e-10

      # EBIC Calculation
      log_lik_term <- n * log(rss / n)
      complexity_penalty <- k * log(n) + 2 * gamma_ebic * lchoose(p, k)
      ebic <- log_lik_term + complexity_penalty

      ebic_scores[k] <- ebic
      rss_scores[k] <- rss
      results_list[[k]] <- fit

      cat(sprintf("Done. RSS=%.2f, EBIC=%.2f\n", rss, ebic))

    } else {
      ebic_scores[k] <- Inf
      rss_scores[k] <- Inf
      results_list[[k]] <- NULL
      cat(sprintf("Failed (Status: %s)\n", fit$status))
    }
  }

  best_k_idx <- which.min(ebic_scores)

  if (length(best_k_idx) == 0 || is.infinite(ebic_scores[best_k_idx])) {
    warning("No valid model found.")
    return(NULL)
  }

  best_model <- results_list[[best_k_idx]]
  best_model$best_k_selected <- best_k_idx
  best_model$all_ebic <- ebic_scores
  best_model$all_rss <- rss_scores

  cat("------------------------------------------------\n")
  cat(sprintf("Winner: k = %d with EBIC = %.2f\n", best_k_idx, ebic_scores[best_k_idx]))

  return(best_model)
}