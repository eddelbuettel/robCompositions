library(testthat)
library(robCompositions)

## cellPcaCoDa is exercised in its valid regime: D >= 3 parts and, for cellwise
## detection, k < D - 1 (a non-empty residual space).  Five-part `expenditures`
## (k = 2 < 4) is used for the structural checks.

test_that("cellPcaCoDa returns correct structure", {
  data(expenditures)
  res <- cellPcaCoDa(expenditures, k = 2)
  expect_s3_class(res, "cellPcaCoDa")
  expect_equal(ncol(res$scores), 2)
  expect_equal(nrow(res$scores), nrow(expenditures))
  expect_equal(nrow(res$loadings), ncol(expenditures))
  expect_equal(ncol(res$loadings), 2)
  expect_equal(dim(res$cellflags), dim(as.matrix(expenditures)))
  expect_true(is.logical(res$cellflags))
  expect_true(is.logical(res$converged))
  expect_true(is.numeric(res$eigenvalues))
  expect_equal(length(res$eigenvalues), 2)
  expect_true(is.numeric(res$center))
  expect_true(is.numeric(res$iterations))
})

test_that("defaults follow Templ (2026): Tukey loss and MCD init", {
  data(expenditures)
  res <- cellPcaCoDa(expenditures, k = 2)
  expect_equal(res$rho, "tukey")
  expect_equal(res$init, "mcd")
})

test_that("eigenvalues_all is returned for variance reporting", {
  data(expenditures)
  res <- cellPcaCoDa(expenditures, k = 2)
  expect_true(is.numeric(res$eigenvalues_all))
  expect_equal(length(res$eigenvalues_all), ncol(expenditures) - 1)
})

test_that("cellPcaCoDa detects contamination in its valid regime (k < D-1)", {
  set.seed(11)
  D <- 8; n <- 300
  z0 <- MASS::mvrnorm(n, mu = rep(0, D - 1), Sigma = diag(D - 1))
  xc <- as.matrix(pivotCoordInv(z0))
  cont <- contaminate_simplex(xc, epsilon = 0.10, delta = 20, seed = 99)
  res <- cellPcaCoDa(cont$x_contaminated, k = 2, init = "classical")
  recall <- mean(res$cellflags[cont$flags])      # flagged | truly contaminated
  fpr    <- mean(res$cellflags[!cont$flags])     # flagged | clean
  expect_gt(recall, fpr)                          # detection carries signal
  expect_lt(fpr, 0.03)                            # false-positive control (Prop.)
})

test_that("clean data produces few flags (false-positive control)", {
  set.seed(123)
  z <- MASS::mvrnorm(200, mu = c(0, 0, 0), Sigma = diag(3))
  x <- pivotCoordInv(z[, 1:3])                     # D = 4, p = 3, k = 2 < p
  res <- cellPcaCoDa(x, k = 2)
  expect_true(mean(res$cellflags) < 0.05)
})

test_that("cellPcaCoDa matches classical PCA on clean data (tight angle)", {
  data(expenditures)
  res_cell <- cellPcaCoDa(expenditures, k = 2)
  res_cla  <- pcaCoDa(expenditures, method = "classical")
  ## with (almost) no flags, cellPcaCoDa reduces to classical ilr PCA
  cos_angle <- abs(det(crossprod(res_cell$loadings[, 1:2],
                                 res_cla$loadings[, 1:2])))
  expect_true(cos_angle > 0.8)
})

test_that("cellPcaCoDa converges", {
  data(expenditures)
  res <- cellPcaCoDa(expenditures, k = 2)
  expect_true(res$converged)
  expect_true(res$iterations < 100)
})

test_that("init = 'classical' runs and converges", {
  data(expenditures)
  res <- cellPcaCoDa(expenditures, k = 2, init = "classical")
  expect_s3_class(res, "cellPcaCoDa")
  expect_equal(res$init, "classical")
  expect_true(res$converged)
})

test_that("cellPcaCoDa falls back to classical init when MCD fails", {
  set.seed(7)
  ## n <= p makes covMcd fail; the estimator must fall back, not error
  z <- MASS::mvrnorm(5, mu = rep(0, 5), Sigma = diag(5))
  x <- pivotCoordInv(z)                            # n = 5, D = 6
  res <- suppressWarnings(cellPcaCoDa(x, k = 2))
  expect_s3_class(res, "cellPcaCoDa")
  expect_equal(res$init, "classical")
})

test_that("cellPcaCoDa works with D=3 (minimum dimension), k=1", {
  set.seed(456)
  z <- MASS::mvrnorm(100, mu = c(0, 0), Sigma = diag(2))
  x <- pivotCoordInv(z)                            # D = 3, p = 2, k = 1 < p
  res <- cellPcaCoDa(x, k = 1)
  expect_s3_class(res, "cellPcaCoDa")
  expect_equal(ncol(res$scores), 1)
  expect_equal(ncol(res$cellflags), 3)
  expect_true(res$converged)
})

test_that("cellPcaCoDa requires at least 3 parts (D >= 3)", {
  x <- constSum(matrix(c(0.3, 0.7, 0.5, 0.5, 0.9, 0.1, 0.4, 0.6, 0.2, 0.8),
                       ncol = 2, byrow = TRUE))
  expect_error(cellPcaCoDa(x, k = 1), "at least 3 parts")
})

test_that("k = D-1 warns that detection is disabled", {
  set.seed(1)
  z <- MASS::mvrnorm(80, mu = c(0, 0), Sigma = diag(2))
  x <- pivotCoordInv(z)                            # D = 3, p = 2
  expect_warning(cellPcaCoDa(x, k = 2), "detection")
})

test_that("cellPcaCoDa works with rho='tukey' and rho='huber'", {
  data(expenditures)
  expect_equal(cellPcaCoDa(expenditures, k = 2, rho = "tukey")$rho, "tukey")
  expect_equal(cellPcaCoDa(expenditures, k = 2, rho = "huber")$rho, "huber")
})

test_that("cellPcaCoDa works with reg > 0 (regularization)", {
  data(expenditures)
  res <- cellPcaCoDa(expenditures, k = 2, reg = 0.1)
  expect_s3_class(res, "cellPcaCoDa")
  expect_true(res$converged)
})

test_that("cellPcaCoDa with maxiter=1 does not converge", {
  data(expenditures)
  res <- cellPcaCoDa(expenditures, k = 2, maxiter = 1)
  expect_false(res$converged)
  expect_equal(res$iterations, 1)
})

test_that("print and summary methods produce output", {
  data(expenditures)
  res <- cellPcaCoDa(expenditures, k = 2)
  expect_output(print(res), "Cellwise-robust PCA")
  expect_output(summary(res), "Importance of components")
})

test_that("contaminate_simplex produces valid compositions", {
  data(arcticLake)
  cont <- contaminate_simplex(arcticLake, epsilon = 0.1, delta = 10, seed = 1)
  rs_orig <- rowSums(as.matrix(arcticLake))
  rs_cont <- rowSums(cont$x_contaminated)
  expect_equal(rs_cont, rs_orig, tolerance = 1e-10)
  expect_true(sum(cont$flags) > 0)
  expect_equal(dim(cont$x_contaminated), dim(as.matrix(arcticLake)))
  expect_equal(dim(cont$flags), dim(as.matrix(arcticLake)))
})

test_that("contaminate_simplex handles epsilon=0", {
  data(arcticLake)
  cont <- contaminate_simplex(arcticLake, epsilon = 0, delta = 10, seed = 1)
  expect_equal(sum(cont$flags), 0)
  expect_equal(cont$x_contaminated, as.matrix(arcticLake))
})

test_that("contaminate_simplex works with type='replacement'", {
  data(arcticLake)
  cont <- contaminate_simplex(arcticLake, epsilon = 0.1, delta = 5,
                              type = "replacement", seed = 10)
  rs_orig <- rowSums(as.matrix(arcticLake))
  rs_cont <- rowSums(cont$x_contaminated)
  expect_equal(rs_cont, rs_orig, tolerance = 1e-10)
})

test_that("contaminate_simplex works with type='additive'", {
  data(arcticLake)
  cont <- contaminate_simplex(arcticLake, epsilon = 0.1, delta = 3,
                              type = "additive", seed = 11)
  rs_orig <- rowSums(as.matrix(arcticLake))
  rs_cont <- rowSums(cont$x_contaminated)
  expect_equal(rs_cont, rs_orig, tolerance = 1e-10)
})

test_that("cellPcaCoDa rejects negative values", {
  x <- matrix(c(1, 2, 3, -1, 5, 6), ncol = 3)
  expect_error(cellPcaCoDa(x), "strictly positive")
})

test_that("cellPcaCoDa rejects NA values", {
  x <- matrix(c(1, 2, 3, NA, 5, 6), ncol = 3)
  expect_error(cellPcaCoDa(x), "missing values")
})

test_that("cellPcaCoDa stops on a limit cycle instead of hitting maxiter (Gjovik)", {
  data(gjovik)
  el <- intersect(c("Al", "Fe", "Ca", "Mg", "Na", "K", "Ti"), colnames(gjovik))
  x <- gjovik[, el]
  x <- x[complete.cases(x), ]
  for (j in seq_len(ncol(x))) {
    z <- x[, j] <= 0
    if (any(z)) x[z, j] <- (2 / 3) * min(x[!z, j])
  }
  x <- constSum(x, const = 100)
  res <- cellPcaCoDa(x, k = 3)
  expect_true(res$converged)        # cycle detection prevents the maxiter hang
  expect_true(res$iterations < 100)
})
