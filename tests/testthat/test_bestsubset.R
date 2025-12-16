library(Matrix)
library(gurobi)
library(testthat)

context("Best Subset Selection")

test_that("best_subset returns correct output", {
  set.seed(123)
  
  # Small synthetic dataset
  X <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
  beta_true <- c(1, 0, -2, 0, 0)
  y <- X %*% beta_true + rnorm(20, sd = 0.5)
  
  k <- 2
  
  # Run best_subset
  fit <- best_subset(X, y, k = k, M = 5, time_limit = 10, verbose = FALSE)
  
  expect_true("beta_vector" %in% names(fit))
  expect_true(length(fit$beta_vector) == ncol(X))
  expect_true(length(fit$selected_indices) <= k)
  expect_true(is.finite(fit$rss))
  expect_true(is.numeric(fit$rss))
  expect_true(is.character(fit$status))
})

test_that("best_subset_auto selects best k and computes MSE", {
  set.seed(123)
  
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
  
  expect_true("beta_vector" %in% names(best_model))
  expect_true(length(best_model$beta_vector) == ncol(X))
  expect_true(best_model$best_k_selected >= 1)
  expect_true(all(is.finite(best_model$all_mse_val)))
  expect_true(all(is.finite(best_model$all_rss_train)))
})