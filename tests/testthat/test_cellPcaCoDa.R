library(testthat)
library(robCompositions)

test_that("cellPcaCoDa returns correct structure", {
  data(arcticLake)
  res <- cellPcaCoDa(arcticLake, k = 2)
  expect_s3_class(res, "cellPcaCoDa")
  expect_equal(ncol(res$scores), 2)
  expect_equal(nrow(res$scores), nrow(arcticLake))
  expect_equal(nrow(res$loadings), ncol(arcticLake))
  expect_equal(ncol(res$loadings), 2)
  expect_equal(dim(res$cellflags), dim(as.matrix(arcticLake)))
  expect_true(is.logical(res$cellflags))
  expect_true(is.logical(res$converged))
  expect_true(is.numeric(res$eigenvalues))
  expect_equal(length(res$eigenvalues), 2)
  expect_true(is.numeric(res$center))
  expect_true(is.numeric(res$iterations))
})

test_that("cellPcaCoDa detects known contamination", {
  data(arcticLake)
  set.seed(42)
  x <- as.matrix(arcticLake)
  ## Contaminate cell [1,2] heavily
  x[1, 2] <- x[1, 2] * 100
  x[1, ] <- x[1, ] / sum(x[1, ]) * sum(arcticLake[1, ])
  res <- cellPcaCoDa(x, k = 2)
  ## Should flag row 1, part 2
  expect_true(res$cellflags[1, 2])
})

test_that("clean data produces no flags (tight threshold)", {
  set.seed(123)
  ## Generate clean logistic-normal data
  z <- MASS::mvrnorm(200, mu = c(0, 0, 0), Sigma = diag(3))
  x <- pivotCoordInv(z[, 1:3])
  res <- cellPcaCoDa(x, k = 2)
  ## False positive rate should be low
  expect_true(mean(res$cellflags) < 0.05)
})

test_that("cellPcaCoDa matches pcaCoDa on clean data (tight angle)", {
  data(expenditures)
  res_cell <- cellPcaCoDa(expenditures, k = 2)
  res_row <- pcaCoDa(expenditures, method = "robust")
  ## Subspace angle between first 2 PCs should be small on clean data
  cos_angle <- abs(det(crossprod(
    res_cell$loadings[, 1:2],
    res_row$loadings[, 1:2]
  )))
  expect_true(cos_angle > 0.8)
})

test_that("cellPcaCoDa converges", {
  data(arcticLake)
  res <- cellPcaCoDa(arcticLake, k = 2)
  expect_true(res$converged)
  expect_true(res$iterations < 100)
})

test_that("cellPcaCoDa works with D=3 (minimum meaningful dimension)", {
  set.seed(456)
  z <- MASS::mvrnorm(100, mu = c(0, 0), Sigma = diag(2))
  x <- pivotCoordInv(z)  # gives D=3 composition
  ## Contaminate one cell
  x[1, 1] <- x[1, 1] * 50
  x[1, ] <- x[1, ] / sum(x[1, ])
  res <- cellPcaCoDa(x, k = 1)
  expect_s3_class(res, "cellPcaCoDa")
  expect_equal(ncol(res$scores), 1)
  expect_true(res$converged)
  ## Should detect contamination in row 1
  expect_true(any(res$cellflags[1, ]))
})

test_that("contaminate_simplex produces valid compositions", {
  data(arcticLake)
  cont <- contaminate_simplex(arcticLake, epsilon = 0.1, delta = 10, seed = 1)
  ## All rows should still sum to the same constant (within tolerance)
  rs_orig <- rowSums(as.matrix(arcticLake))
  rs_cont <- rowSums(cont$x_contaminated)
  expect_equal(rs_cont, rs_orig, tolerance = 1e-10)
  ## Some cells should be flagged
  expect_true(sum(cont$flags) > 0)
  ## Output dimensions match
  expect_equal(dim(cont$x_contaminated), dim(as.matrix(arcticLake)))
  expect_equal(dim(cont$flags), dim(as.matrix(arcticLake)))
})

test_that("contaminate_simplex handles epsilon=0", {
  data(arcticLake)
  cont <- contaminate_simplex(arcticLake, epsilon = 0, delta = 10, seed = 1)
  ## No flags
  expect_equal(sum(cont$flags), 0)
  ## Data unchanged
  expect_equal(cont$x_contaminated, as.matrix(arcticLake))
})
