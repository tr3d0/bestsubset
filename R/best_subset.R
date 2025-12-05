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
#' @importFrom Matrix Diagonal Matrix
#' @importFrom gurobi gurobi
#' @export
best_subset <- function(X, y, k, M = 15, 
                        time_limit = 300, 
                        mip_gap = 0.01,
                        verbose = TRUE,
                        threads = 0,
                        mip_focus = 3,
                        presolve = 2,
                        ...) {
  if (!requireNamespace("gurobi", quietly = TRUE)) {
    stop("Package 'gurobi' is required but not installed. Install it as described in README.")
  }
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("`X` must be numeric matrix.")
  y <- as.numeric(y)
  if (nrow(X) != length(y)) stop("`X` and `y` have incompatible dimensions.")
  if (!is.numeric(k) || k < 0 || k > ncol(X)) stop("`k` must be integer between 0 and p.")

  # 1. Data Preparation (Standardization)
  X_scaled <- scale(X, center = TRUE, scale = TRUE)
  y_scaled <- scale(y, center = TRUE, scale = FALSE)
  
  scale_X <- attr(X_scaled, "scaled:scale")
  scale_X[scale_X == 0] <- 1  # evita divisione per zero (coefficienti resteranno 0)
  center_X <- attr(X_scaled, "scaled:center")
  mean_y <- mean(y)
  
  n <- nrow(X_scaled)
  p <- ncol(X_scaled)
  
  # 2. Quadratic Objective
  XtX <- crossprod(X_scaled)
  Xty <- crossprod(X_scaled, y_scaled)
  XtX <- as(XtX, "generalMatrix")
  
  # Q matrix: [2X'X, 0; 0, 0]
  # Gurobi's quadratic term is 1/2 * x' Q x, so set Q = 2 * X'X to get x'X'X x
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
  A3 <- Matrix::Matrix(c(rep(0, p), rep(1, p)), nrow = 1, sparse = TRUE)
    
  A <- rbind(A1, A2, A3)
  A <- as(A, "generalMatrix")
  rhs <- c(rep(0, 2 * p), k)
  sense <- c(rep("<=", 2 * p), "=")
  
  # 4. Gurobi Model
  model <- list()
  model$Q <- as(Q, "generalMatrix")
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
    if (!is.null(colnames(X))) {
      names(beta_final) <- colnames(X)
    } else {
      names(beta_final) <- paste0("V", seq_len(p))
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

#' Best Subset Selection with Validation Set
#'
#' @description
#' Runs the best subset selection algorithm for a range of subset sizes (from 1 to \code{k_max})
#' and selects the globally optimal model minimizing the Mean Square Error (MSE).
#'
#' @details
#' This function acts as a wrapper around \code{\link{best_subset}}. It iterates through
#' all possible subset sizes up to \code{k_max}. For each size \code{k}, it fits the best subset
#' model using the MIQP solver.
#' Runs the best subset selection algorithm for a range of subset sizes using MIQP.
#' It fits the models on the training set (X, y) and selects the best k by minimizing 
#' the Mean Squared Error (MSE) on the provided validation/test set (X_val, y_val).
#'
#' @param X Numeric matrix. The design matrix (predictors) - TRAINING.
#' @param y Numeric vector. The response vector - TRAINING.
#' @param X_val Numeric matrix. The design matrix (predictors) - VALIDATION/TEST.
#' @param y_val Numeric vector. The response vector - VALIDATION/TEST.
#' @param k_max Integer. Maximum subset size to consider. Default is 10.
#' @param time_limit_per_k Numeric. Time limit per k in seconds. Default is 100.
#' @param mip_gap Numeric. The relative MIP optimality gap. Default is 0.01.
#' @param plot Logical. If TRUE, generates an Elbow plot (RSS) and an MSE plot. Default is TRUE.
#' @param verbose Logical. If TRUE, shows progress. Default is FALSE.
#' @param ... Additional parameters passed to \code{\link{best_subset}} (e.g., \code{threads}).
#'
#' @export
best_subset_auto <- function(X, y, X_val, y_val, 
                             k_max = 10,
                             time_limit_per_k = 100, 
                             mip_gap = 0.01, 
                             plot = TRUE,
                             verbose = FALSE, 
                             ...) {
  
  n <- nrow(X)
  p <- ncol(X)
  if (k_max > p) k_max <- p
  
  results_list <- list()
  mse_scores <- rep(Inf, k_max)
  rss_scores <- rep(Inf, k_max)
  
  message("Starting Best Subset Selection (k_max = ", k_max, ")...")
  message("------------------------------------------------")
  
  for (k in 1:k_max) {
    message(sprintf("Solving for k = %2d ... ", k), appendLF = FALSE)
    
    fit <- best_subset(X, y, k = k, 
                       time_limit = time_limit_per_k, 
                       verbose = FALSE, # Suppress inner Gurobi output to keep console clean
                       ...)
    
    if (fit$status %in% c("OPTIMAL", "TIME_LIMIT", "SUBOPTIMAL")) {
      
      rss_scores[k] <- fit$rss
      
      
      intercept <- fit$coefficients["Intercept"]
      
      y_pred_val <- intercept + (X_val %*% fit$beta_vector)
      
      mse_val <- mean((y_val - y_pred_val)^2)
      
      mse_scores[k] <- mse_val
      
      results_list[[k]] <- fit
      
      message("Done. Val MSE=", sprintf("%.4f", mse_val))
      
    } else {
      message("Failed (Status: ", fit$status, ")")
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
  
  message("------------------------------------------------")
  message("Winner: k = ", best_k_idx, " with MSE = ", sprintf("%.4f", mse_scores[best_k_idx]))

  
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