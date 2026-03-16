library(testthat)
library(robCompositions)

test_that("closure propagation matches theoretical prediction", {
  ## Test that contaminating part j and re-closing produces
  ## the predicted perturbation in clr coordinates
  D <- 5
  x <- c(0.2, 0.3, 0.1, 0.15, 0.25)  # composition on S^5
  delta <- 3.0
  j <- 2  # contaminate part 2

  ## Contaminate and re-close
  x_cont <- x
  x_cont[j] <- x_cont[j] * delta
  x_cont <- x_cont / sum(x_cont)

  ## Observed clr perturbation
  clr_clean <- log(x) - mean(log(x))
  clr_cont <- log(x_cont) - mean(log(x_cont))
  delta_clr_obs <- clr_cont - clr_clean

  ## Theoretical prediction (Theorem 3.1 from propagation_theorem.tex)
  delta_clr_theory <- rep(-log(delta) / D, D)
  delta_clr_theory[j] <- (D - 1) / D * log(delta)

  expect_equal(delta_clr_obs, delta_clr_theory, tolerance = 1e-12)
})

test_that("ilr perturbation is rank-1", {
  D <- 10
  V <- orthbasis(D)$V
  j <- 3
  delta <- 5.0

  ## Theoretical clr perturbation for contamination of part j
  delta_clr <- rep(-log(delta) / D, D)
  delta_clr[j] <- (D - 1) / D * log(delta)

  ## Map to ilr
  delta_ilr <- as.numeric(t(V) %*% delta_clr)

  ## This should be proportional to V^T(e_j - 1/D * 1)
  expected_direction <- as.numeric(t(V) %*% (diag(D)[, j] - rep(1/D, D)))
  expected_direction <- expected_direction / sqrt(sum(expected_direction^2))
  delta_ilr_norm <- delta_ilr / sqrt(sum(delta_ilr^2))

  ## Inner product should be exactly 1 (parallel vectors)
  expect_equal(abs(sum(delta_ilr_norm * expected_direction)),
               1.0, tolerance = 1e-12)
})

test_that("propagation works for multiple contaminated parts", {
  D <- 6
  V <- orthbasis(D)$V
  x <- rep(1/D, D)  # uniform composition
  delta_j <- 3.0
  delta_k <- 2.0
  j <- 2; k <- 5

  ## Contaminate two parts simultaneously
  x_cont <- x
  x_cont[j] <- x_cont[j] * delta_j
  x_cont[k] <- x_cont[k] * delta_k
  x_cont <- x_cont / sum(x_cont)

  ## Observed ilr perturbation
  z_clean <- as.numeric(t(V) %*% (log(x) - mean(log(x))))
  z_cont <- as.numeric(t(V) %*% (log(x_cont) - mean(log(x_cont))))
  delta_ilr_obs <- z_cont - z_clean

  ## Theoretical: sum of two rank-1 perturbations
  d_j <- as.numeric(t(V) %*% (diag(D)[, j] - 1/D))
  d_k <- as.numeric(t(V) %*% (diag(D)[, k] - 1/D))
  delta_ilr_theory <- log(delta_j) * d_j + log(delta_k) * d_k

  expect_equal(delta_ilr_obs, delta_ilr_theory, tolerance = 1e-12)
})

test_that("single-part contamination affects all ilr coordinates", {
  D <- 8
  V <- orthbasis(D)$V

  for (j in 1:D) {
    d_j <- as.numeric(t(V) %*% (diag(D)[, j] - 1/D))
    ## All D-1 ilr coordinates should be non-zero
    expect_true(all(abs(d_j) > 1e-14),
                info = sprintf("Part %d: some ilr coordinates are zero", j))
  }
})
