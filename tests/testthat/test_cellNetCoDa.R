## cellNetCoDa: end-to-end estimator (skipped on CRAN -- StARS subsampling is
## slow) plus fast S3-method coverage on a constructed object.

test_that("cellNetCoDa runs end-to-end for both back ends and returns a graph", {
  skip_on_cran()
  skip_if_not_installed("glmnet")
  skip_if_not_installed("glasso")
  set.seed(9)
  D <- 8; n <- 150
  X <- matrix(rpois(n * D, lambda = 6), n, D); X[3, 2] <- 0   # counts with a zero
  res_mb <- cellNetCoDa(X, method = "mb",
                        lambda_grid = c(0.4, 0.2, 0.1, 0.05), B = 8, k = 2)
  res_gl <- cellNetCoDa(X, method = "glasso",
                        rho_grid = c(0.4, 0.2, 0.1, 0.05), B = 8, k = 2)
  for (res in list(res_mb, res_gl)) {
    expect_s3_class(res, "cellNetCoDa")
    expect_equal(dim(res$adj), c(D, D))
    expect_equal(res$adj, t(res$adj))             # symmetric
    expect_true(all(diag(res$adj) == 0))          # zero diagonal
    expect_true(is.matrix(res$weights) && all(res$weights >= 0 & res$weights <= 1))
  }
  ## S3 methods on a real fit
  expect_output(print(res_mb), "Cellwise-robust network")
  sm <- summary(res_mb)
  expect_true(is.list(sm) && all(c("degree", "flagrate", "edges") %in% names(sm)))
})

test_that("cellNetCoDa S3 methods work on a constructed object (CRAN-fast)", {
  set.seed(1)
  D <- 6; n <- 30
  A <- matrix(0, D, D)
  A[upper.tri(A)] <- rbinom(D * (D - 1) / 2, 1, 0.35)
  A <- A + t(A); diag(A) <- 0
  dimnames(A) <- list(paste0("P", 1:D), paste0("P", 1:D))
  flags <- matrix(runif(n * D) < 0.08, n, D)
  obj <- structure(list(adj = A, taxa = paste0("P", 1:D),
                        n_edges = sum(A[upper.tri(A)]),
                        weights = ifelse(flags, 0.05, 1),
                        penalty = 0.1, instability = c(0.002, 0.02, 0.04, 0.09),
                        grid = c(0.4, 0.2, 0.1, 0.05), thresh = 0.05,
                        method = "mb", flags = flags),
                   class = "cellNetCoDa")

  expect_output(print(obj), "density")
  sm <- summary(obj)
  expect_equal(length(sm$degree), D)
  expect_s3_class(sm$edges, "data.frame")

  pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(obj)                 # which = 1 (network)
  plot(obj, which = 2)      # contamination heatmap
  plot(obj, which = 3)      # StARS instability path
  expect_error(plot(obj, which = 9), "must be 1")
})
