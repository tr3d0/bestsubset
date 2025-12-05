# bestsubset: Exact Best Subset Selection via MIQP

**bestsubset** is an R package that implements an exact algorithm for **Best Subset Selection** in linear regression. It formulates the variable selection problem as a Mixed-Integer Quadratic Programming (MIQP) problem and solves it using the **Gurobi** optimizer.

Unlike stepwise methods (forward/backward) or Lasso (L1 regularization), this package solves the $L_0$ regularization problem **exactly**, finding the globally optimal subset of predictors for a given size $k$.

## Features

- **Exact Global Optimum**: Uses MIQP to find the absolute best subset of size $k$ (no heuristics).
- **Big-M Constraints**: Implements tight Big-M constraints on regression coefficients.
- **Auto-Tuning (Validation set)**: Automatically selects the optimal model size ($k$) using the **Validation Set approach** (minimizing MSE).
- **Robust Wrapper**: Built-in scaling, time limits, and optimality gap controls.

## Prerequisites

This package requires the **Gurobi Optimizer** software and its R interface.

1.  **Install Gurobi Software**: Download and install Gurobi from [gurobi.com](https://www.gurobi.com/). You will need a valid license (academic licenses are free).
2.  **Install the Gurobi R Package**:
    ```r
    install.packages("gurobi", repos = "[https://gurobi.github.io/drat](https://gurobi.github.io/drat)")
    # OR follow the instructions in your Gurobi installation folder
    ```

## Installation

Once Gurobi is ready, you can install **bestsubset** directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("tr3d0/bestsubset")
```

## Usage

### 1. Simulated Data
Let's create a synthetic dataset with 100 observations and 20 predictors, where only the first 3 are truly significant.

```r
library(bestsubset)

set.seed(123)
n <- 100
p <- 20
X <- matrix(rnorm(n * p), n, p)
beta_true <- c(3, 2, 1.5, rep(0, p - 3)) # Only first 3 are non-zero
y <- X %*% beta_true + rnorm(n)
```

### 2. Automatic Selection (Recommended)
Use `best_subset_auto()` to check multiple subset sizes and pick the winner via Validation MSE. You need to provide a training and a validation set.

```r
# Split data into Train (70%) and Validation (30%)
train_idx <- sample(1:n, size = 0.7 * n)
X_train <- X[train_idx, ]; y_train <- y[train_idx]
X_val   <- X[-train_idx, ]; y_val   <- y[-train_idx]

# Search for best subsets from k=1 to k=10
# Time limit: 100s per k
fit <- best_subset_auto(X_train, y_train, X_val, y_val, 
                        k_max = 10, 
                        verbose = FALSE)

# Print the winner
print(fit$coefficients)

# Check which variables were selected (Should be 1, 2, 3)
print(fit$selected_indices)

# View the best k
cat("Optimal subset size:", fit$best_k_selected)
```

### 3. Fixed Size Selection
If you know exactly how many variables you want (e.g., $k=3$) or want to use the full dataset for a specific fit:

```r
# Force selection of exactly 3 variables
fit_k3 <- best_subset(X, y, k = 3, time_limit = 60)

print(fit_k3$coefficients)
```

## Parameters

- **`time_limit`**: Maximum time (in seconds) allowed for the solver.
- **`mip_gap`**: Stop optimization when the solution is within this % of the theoretical optimum (default 1% or `0.01`).
- **`verbose`**: Set to `TRUE` to see the live Gurobi solver logs.

## License

MIT License. See `LICENSE` file for details.

## Author

**Edoardo Lanzetti**