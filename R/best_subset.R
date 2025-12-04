#' Best Subset Selection using Gurobi
#'
#' Solves the best subset selection problem using Mixed-Integer Quadratic Programming (MIQP).
#'
#' @param X Numeric matrix. The design matrix (predictors).
#' @param y Numeric vector. The response vector.
#' @param k Integer. The exact number of variables to select.
#' @param M Numeric. Big-M parameter for constraints. Default is 100.
#' @param time_limit Numeric. Time limit for the solver in seconds. Default is 60.
#'
#' @return A list containing:
#' \item{coefficients}{A named vector of estimated coefficients (including intercept).}
#' \item{selected_indices}{Indices of the selected variables.}
#' \item{rss}{The Residual Sum of Squares of the fit.}
#' \item{status}{The status of the Gurobi optimization (e.g., "OPTIMAL").}
#'
#' @importFrom Matrix crossprod Diagonal Matrix as
#' @importFrom stats start
#' @import gurobi
#' @export
best_subset <- function(X, y, k, M = 100, time_limit = 60) {

  # --- 1. Data Prep (Centering) ---
  X_scaled <- scale(X, center = TRUE, scale = FALSE)
  y_scaled <- scale(y, center = TRUE, scale = FALSE)

  mean_y <- mean(y)
  mean_X <- colMeans(X)

  n <- nrow(X_scaled)
  p <- ncol(X_scaled)

  # --- 2. Quadratic Objective ---
  XtX <- crossprod(X_scaled)
  Xty <- crossprod(X_scaled, y_scaled)

  # Matrix package compatibility fix
  XtX <- as(XtX, "generalMatrix")

  # Matrix Q: [2X'X, 0; 0, 0]
  Q <- Matrix(0, 2 * p, 2 * p, sparse = TRUE)
  Q[1:p, 1:p] <- 2 * XtX

  # Linear vector: -2X'y for beta, 0 for z
  obj <- c(-2 * Xty, rep(0, p))

  # --- 3. Constraints ---
  # Big-M: |beta| <= M*z
  A1 <- cbind(Diagonal(p), Diagonal(p, x = -M))          # beta - Mz <= 0
  A2 <- cbind(Diagonal(p, x = -1), Diagonal(p, x = -M))  # -beta - Mz <= 0

  # Cardinality: sum(z) = k
  A3 <- c(rep(0, p), rep(1, p))

  A <- rbind(A1, A2, A3)
  rhs <- c(rep(0, 2 * p), k)
  sense <- c(rep("<", 2 * p), "=")

  # --- 4. Gurobi Model ---
  model <- list()
  model$Q <- Q
  model$obj <- obj
  model$A <- A
  model$rhs <- rhs
  model$sense <- sense
  model$vtype <- c(rep("C", p), rep("B", p))

  params <- list(OutputFlag = 1, TimeLimit = time_limit)

  # --- 5. Solve and Output ---
  result <- gurobi(model, params)

  if (result$status %in% c("OPTIMAL", "TIME_LIMIT")) {
    beta_hat <- result$x[1:p]
    z_hat <- result$x[(p + 1):(2 * p)]

    intercept <- mean_y - sum(mean_X * beta_hat)

    return(list(
      coefficients = c(Intercept = intercept, beta_hat),
      selected_indices = which(z_hat > 0.5),
      rss = result$objval + sum(y_scaled^2),
      status = result$status
    ))
  } else {
    return(list(status = result$status))
  }
}