#' Best Subset Selection via MIQP (Fixed k)
#'
#' @description
#' Solves the best subset selection problem for a fixed subset size \code{k} using
#' Mixed-Integer Quadratic Programming (MIQP).
#'
#' @param X Numeric matrix. The design matrix (predictors).
#' @param y Numeric vector. The response vector.
#' @param k Integer. The exact number of variables to select.
#' @param M Numeric. Big-M parameter for constraints. Default is 15.
#' @param time_limit Numeric. Time limit for the solver in seconds. Default is 300.
#' @param mip_gap Numeric. The relative MIP optimality gap. Default is 0.01 (1%).
#' @param verbose Logical. If TRUE, shows Gurobi console output. Default is TRUE.
#' @param threads Integer. Number of threads to use (0 = all available). Default is 0.
#' @param ... Additional parameters passed directly to the Gurobi solver.
#'
#' @return A list containing coefficients, RSS, and status.
#' @importFrom Matrix crossprod Diagonal Matrix as
#' @import gurobi
#' @export
best_subset <- function(X, y, k, M = 15,
                        time_limit = 300,
                        mip_gap = 0.01,
                        verbose = TRUE,
                        threads = 0,
                        ...) {

  X_scaled <- scale(X, center = TRUE, scale = TRUE)
  y_scaled <- scale(y, center = TRUE, scale = FALSE)

  scale_X <- attr(X_scaled, "scaled:scale")
  center_X <- attr(X_scaled, "scaled:center")
  mean_y <- mean(y)

  n <- nrow(X_scaled)
  p <- ncol(X_scaled)
  
  XtX <- crossprod(X_scaled)
  Xty <- crossprod(X_scaled, y_scaled)
  XtX <- as(XtX, "generalMatrix")
  
  Q <- Matrix(0, 2 * p, 2 * p, sparse = TRUE)
  Q[1:p, 1:p] <- 2 * XtX
  
  obj <- c(-2 * Xty, rep(0, p))
  
  A1 <- cbind(Diagonal(p), Diagonal(p, x = -M))
  A2 <- cbind(Diagonal(p, x = -1), Diagonal(p, x = -M))
  A3 <- c(rep(0, p), rep(1, p))
  
  A <- rbind(A1, A2, A3)
  rhs <- c(rep(0, 2 * p), k)
  sense <- c(rep("<", 2 * p), "=")
  
  model <- list()
  model$Q <- Q
  model$obj <- obj
  model$A <- A
  model$rhs <- rhs
  model$sense <- sense
  model$vtype <- c(rep("C", p), rep("B", p))
  
  params <- list(
    OutputFlag = as.numeric(verbose),
    TimeLimit = time_limit,
    MIPGap = mip_gap,
    Threads = threads
  )
  
  extra_params <- list(...)
  if (length(extra_params) > 0) params <- c(params, extra_params)
  
  result <- gurobi(model, params)
  
  if (result$status %in% c("OPTIMAL", "TIME_LIMIT", "SUBOPTIMAL")) {
    beta_hat_scaled <- result$x[1:p]
    z_hat <- result$x[(p + 1):(2 * p)]
    selected_vars <- which(z_hat > 0.5)
    
    beta_final <- rep(0, p)
    if(length(selected_vars) > 0) {
      beta_final[selected_vars] <- beta_hat_scaled[selected_vars] / scale_X[selected_vars]
    }
    
    intercept <- mean_y - sum(center_X * beta_final)
    
    y_pred <- intercept + X %*% beta_final
    rss_real <- sum((y - y_pred)^2)
    
    return(list(
      coefficients = c(Intercept = intercept, beta_final),
      beta_vector = beta_final,
      intercept = intercept,
      selected_indices = selected_vars,
      rss = rss_real,
      gap = result$mipgap,
      status = result$status
    ))
  } else {
    return(list(status = result$status, rss = Inf, selected_indices = integer(0)))
  }
}

#' Best Subset Selection with Validation Set (MSE Selection)
#'
#' @description
#' Runs best subset selection for k=1 to k_max on a Training set, 
#' and selects the optimal k by minimizing Mean Squared Error (MSE) on a Validation set.
#'
#' @param X_train Numeric matrix. Training predictors.
#' @param y_train Numeric vector. Training response.
#' @param X_val Numeric matrix. Validation predictors.
#' @param y_val Numeric vector. Validation response.
#' @param k_max Integer. Maximum subset size to consider.
#' @param time_limit_per_k Numeric. Time limit per k in seconds.
#' @param plot Logical. If TRUE, plots RSS (Train) and MSE (Validation).
#' @param verbose Logical. If TRUE, shows loop progress.
#' @param ... Additional parameters passed to \code{\link{best_subset}}.
#'
#' @return A list containing the best model and history of metrics.
#' @importFrom graphics par plot points abline legend grid lines
#' @export
best_subset_auto <- function(X_train, y_train, X_val, y_val, 
                                        k_max = 10,
                                        time_limit_per_k = 100, 
                                        plot = TRUE,
                                        verbose = TRUE, 
                                        ...) {
  
  n_val <- length(y_val)
  p <- ncol(X_train)
  if(k_max > p) k_max <- p
  
  results_list <- list()
  mse_scores <- rep(Inf, k_max)
  rss_scores <- rep(Inf, k_max)
  
  if(verbose) {
    cat("Starting Best Subset Selection (Validation Set Approach)\n")
    cat(sprintf("Train n=%d, Val n=%d, p=%d, k_max=%d\n", nrow(X_train), n_val, p, k_max))
    cat("------------------------------------------------\n")
  }
  
  for (k in 1:k_max) {
    if(verbose) cat(sprintf("  k = %2d ... ", k))
    
    fit <- best_subset(X_train, y_train, k = k, 
                       time_limit = time_limit_per_k, 
                       verbose = FALSE, 
                       ...)
    
    if (fit$status %in% c("OPTIMAL", "TIME_LIMIT", "SUBOPTIMAL")) {
      
      rss_scores[k] <- fit$rss
      pred_val <- fit$intercept + (X_val %*% fit$beta_vector)
      
      mse_val <- mean((y_val - pred_val)^2)
      mse_scores[k] <- mse_val
      
      results_list[[k]] <- fit
      
      if(verbose) cat(sprintf("Done. Train RSS=%.1f, Val MSE=%.4f\n", fit$rss, mse_val))
      
    } else {
      if(verbose) cat(sprintf("Failed (Status: %s)\n", fit$status))
    }
  }
  
  best_k_idx <- which.min(mse_scores)
  
  if (length(best_k_idx) == 0 || is.infinite(mse_scores[best_k_idx])) {
    warning("No valid model found.")
    return(NULL)
  }
  
  best_model <- results_list[[best_k_idx]]
  best_model$best_k_selected <- best_k_idx
  best_model$all_mse_val <- mse_scores
  best_model$all_rss_train <- rss_scores
  
  if(verbose) {
    cat("------------------------------------------------\n")
    cat(sprintf("Winner: k = %d with MSE = %.4f\n", best_k_idx, mse_scores[best_k_idx]))
  }
  
  if (plot) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    
    plot(1:k_max, mse_scores, type = "b", pch = 19, col = "darkred",
         xlab = "Subset Size (k)", ylab = "Mean Squared Error (MSE)",
         main = "Validation Error", xaxt = "n")
    axis(1, at = 1:k_max)
    grid()
    
    points(best_k_idx, mse_scores[best_k_idx], col = "blue", pch = 8, cex = 2, lwd = 2)
    abline(v = best_k_idx, col = "blue", lty = 2)
    legend("topright", legend = c("Best k (Min MSE)"), 
           col = "blue", pch = 8, lty = 2, bty = "n", cex = 0.8)
  }
  
  return(best_model)
}