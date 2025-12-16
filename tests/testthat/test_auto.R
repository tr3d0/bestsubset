library(testthat)
library(Matrix)
library(gurobi)

context("Best Subset Auto Selection")

test_that("best_subset_auto returns valid model and computes MSE", {
  set.seed(123)
  
  # Small synthetic dataset
  X <- matrix(rnorm(30 * 5), nrow = 30, ncol = 5)
  beta_true <- c(1, 0, -2, 0, 0)
  y <- X %*% beta_true + rnorm(30, sd = 0.5)
  
  # Split train/validation
  train_idx <- 1:20
  val_idx <- 21:30
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_val <- X[val_idx, ]
  y_val <- y[val_idx]
  
  # Run best_subset_auto
  best_model <- best_subset_auto(X_train, y_train, X_val, y_val, 
                                 k_max = 3, time_limit_per_k = 10, verbose = FALSE)
  
  # Basic checks
  expect_true("beta_vector" %in% names(best_model))
  expect_true(length(best_model$beta_vector) == ncol(X))
  expect_true(best_model$best_k_selected >= 1 && best_model$best_k_selected <= 3)
  expect_true(all(is.finite(best_model$all_mse_val)))
  expect_true(all(is.finite(best_model$all_rss_train)))
  
  # Optional: check that at least one variable was selected
  expect_true(length(best_model$selected_indices) > 0)
})