library(testthat)
library(bestsubset)
library(Matrix)

context("Best Subset Selection Tests")

test_that("best_subset recovers true support in noiseless setting", {
  # Salta il test se Gurobi non è installato (fondamentale per GitHub Actions)
  skip_if_not_installed("gurobi")
  
  set.seed(123)
  n <- 50
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  # Solo le prime 3 variabili sono attive e ben separate
  beta_true <- c(5, -3, 2, rep(0, p - 3)) 
  y <- X %*% beta_true # Nessun rumore
  
  # Chiediamo esattamente k=3
  fit <- best_subset(X, y, k = 3, time_limit = 20, verbose = FALSE)
  
  expect_equal(fit$status, "OPTIMAL")
  expect_setequal(fit$selected_indices, c(1, 2, 3))
  expect_lt(fit$rss, 1e-4) # RSS deve essere quasi zero
})

test_that("best_subset_auto runs and selects a valid model", {
  skip_if_not_installed("gurobi")
  
  set.seed(42)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- c(3, 3, rep(0, p - 2))
  y <- X %*% beta_true + rnorm(n, sd = 0.5)
  
  # Split manuale
  train <- 1:80
  X_train <- X[train, ]
  y_train <- y[train]
  X_val <- X[-train, ]
  y_val <- y[-train]
  
  # Testiamo il wrapper automatico (che usa il warm start internamente)
  fit_auto <- best_subset_auto(X_train, y_train, X_val, y_val, 
                               k_max = 5, 
                               time_limit_per_k = 5,
                               plot = FALSE, 
                               verbose = FALSE)
  
  expect_type(fit_auto, "list")
  expect_true("best_k_selected" %in% names(fit_auto))
  # Ci aspettiamo che trovi un modello sensato (k tra 1 e 5)
  expect_true(fit_auto$best_k_selected >= 1 && fit_auto$best_k_selected <= 5)
  expect_true(length(fit_auto$beta_vector) == p)
})

test_that("Warm start arguments are accepted without error", {
  skip_if_not_installed("gurobi")
  
  n <- 30; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  
  # Primo fit per avere dei valori
  fit1 <- best_subset(X, y, k = 1, time_limit = 5, verbose = FALSE)
  
  # Secondo fit usando i valori del primo come warm start
  # Se il codice fosse rotto, questo darebbe errore
  expect_error(
    best_subset(X, y, k = 2, 
                start_beta = fit1$beta_scaled, 
                start_z = fit1$z_vector,
                time_limit = 5, verbose = FALSE),
    NA # NA significa "ci aspettiamo NO error"
  )
})

test_that("Input validation works correctly", {
  skip_if_not_installed("gurobi")
  
  n <- 20; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  
  expect_error(best_subset(X, y, k = -1), "must be integer")
  expect_error(best_subset(X, y, k = p + 1), "must be integer")
  expect_error(best_subset(X, y[-1], k = 2), "incompatible dimensions")
})

test_that("Big-M warning triggers correctly", {
  skip_if_not_installed("gurobi")
  
  set.seed(1)
  n <- 30; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  # Creiamo un y che richiede un coefficiente GIGANTE (> M=15)
  y <- 1000 * X[,1] + rnorm(n) 
  
  # Ci aspettiamo un warning
  expect_warning(best_subset(X, y, k = 1, M = 15, verbose = FALSE), 
                 "Coefficient hit Big-M bound")
})