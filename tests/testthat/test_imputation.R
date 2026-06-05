## Characterization + regression tests for the detection-limit / rounded-zero
## imputation family: imputeBDL (canonical), impRZalr, impAll, and the
## deprecated aliases imputeBDLs / impRZilr.
##
## Golden values were captured from the current implementation on deterministic
## paths (method = "lm", no bootstrap) so the imputeBDL rename and the cleanup
## can be proven behaviour-preserving. See plan: imputeBDLs review (2026-06-04).

context("rounded-zero / detection-limit imputation")

test_that("imputeBDL (arcticLake, lm) reproduces reference imputation", {
  data(arcticLake)
  x <- arcticLake
  z_silt <- x[, 2] < 44          # silt below its detection limit
  z_sand <- x[, 1] < 5           # sand below its detection limit (dl = 5)
  x[z_silt, 2] <- 0
  r <- imputeBDL(x, dl = c(5, 44, 0), eps = 0.01, method = "lm", variation = FALSE)

  expect_s3_class(r, "replaced")
  expect_equal(dim(r$x), c(39L, 3L))
  expect_equal(r$method, "lm")
  ## golden values (imputed silt, first four rows)
  expect_equal(round(unname(r$x[1:4, "silt"]), 5),
               c(23.84795, 24.10034, 36.08242, 30.54653),
               tolerance = 1e-4)
  ## invariant: every imputed cell stays at or below its detection limit
  expect_true(all(r$x[z_silt, "silt"] <= 44 + 1e-6))
  expect_true(all(r$x[z_sand, "sand"] <= 5 + 1e-6))
  ## invariant: observed (non-zero) cells are left untouched
  expect_equal(r$x[!z_silt, "silt"], arcticLake[!z_silt, "silt"], tolerance = 1e-8)
  ## wind (zero index) is aligned with the returned original-order columns
  expect_equal(colnames(r$wind), colnames(r$x))
})

test_that("imputeBDL (simulated, lm) reproduces reference imputation", {
  set.seed(1)
  p <- 10; n <- 50; k <- 2
  Tm <- matrix(rnorm(n * k), ncol = k)
  Bm <- matrix(runif(p * k, -1, 1), ncol = k)
  XE <- Tm %*% t(Bm) + matrix(rnorm(n * p, 0, 0.1), ncol = p)
  d <- data.frame(pivotCoordInv(XE))
  pc <- ncol(d)                                   # = 11 parts
  DL <- numeric(pc)
  for (j in seq(1, pc, 2)) DL[j] <- quantile(d[, j], 0.06)
  for (j in seq_len(pc)) d[d[, j] < DL[j], j] <- 0

  r <- imputeBDL(d, dl = DL, eps = 1, method = "lm", variation = FALSE, R = 10)

  expect_s3_class(r, "replaced")
  expect_equal(dim(r$x), c(50L, 11L))
  expect_equal(sum(r$wind), 18L)
  expect_equal(sum(as.matrix(r$x)), 49.979132, tolerance = 1e-3)
  expect_equal(round(r$x[1, 1], 6), 0.143137, tolerance = 1e-5)
  for (j in which(DL > 0)) {
    zj <- which(d[, j] == 0)
    if (length(zj)) expect_true(all(r$x[zj, j] <= DL[j] + 1e-8))
  }
})

test_that("impRZalr (arcticLake, pos = 3) reproduces reference imputation", {
  data(arcticLake)
  x <- arcticLake
  x[x[, 1] < 5, 1] <- 0
  x[x[, 2] < 47, 2] <- 0
  r <- impRZalr(x, pos = 3, dl = c(5, 47), eps = 0.05)

  expect_true(is.list(r))
  expect_true(all(c("x", "iter") %in% names(r)))
  expect_equal(dim(r$x), c(39L, 3L))
  expect_equal(round(unname(r$x[1:2, "silt"]), 5),
               c(22.71156, 22.98108),
               tolerance = 1e-4)
})

test_that("impAll returns imputed data, not NULL, for rounded zeros (P0-1 regression)", {
  skip_on_cran()
  set.seed(1)
  p <- 10; n <- 50; k <- 2
  Tm <- matrix(rnorm(n * k), ncol = k)
  Bm <- matrix(runif(p * k, -1, 1), ncol = k)
  d <- data.frame(pivotCoordInv(Tm %*% t(Bm) + matrix(rnorm(n * p, 0, 0.1), ncol = p)))
  ## impAll convention: values below the detection limit are coded as -dl
  thr <- as.numeric(quantile(d[, 1], 0.10))
  d[d[, 1] < thr, 1] <- -thr

  res <- impAll(d)

  expect_false(is.null(res))                      # was NULL before the res$x fix
  expect_equal(dim(res), c(n, p + 1L))
  expect_true(all(res > 0))                        # negatives/zeros replaced
})

test_that("deprecated imputeBDLs / impRZilr warn and match imputeBDL", {
  data(arcticLake)
  x <- arcticLake
  x[x[, 2] < 44, 2] <- 0

  expect_warning(impRZilr(x, dl = c(5, 44, 0), eps = 0.01, method = "lm"),
                 "imputeBDL")   # locale-independent (message names the replacement)

  ref  <- imputeBDL(x, dl = c(5, 44, 0), eps = 0.01, method = "lm", variation = FALSE)
  wrap <- suppressWarnings(impRZilr(x, dl = c(5, 44, 0), eps = 0.01, method = "lm"))
  expect_equal(wrap$x, ref$x, tolerance = 1e-10)   # wrapper == canonical

  set.seed(1)
  p <- 10; n <- 50; k <- 2
  Tm <- matrix(rnorm(n * k), ncol = k)
  Bm <- matrix(runif(p * k, -1, 1), ncol = k)
  d <- data.frame(pivotCoordInv(Tm %*% t(Bm) + matrix(rnorm(n * p, 0, 0.1), ncol = p)))
  pc <- ncol(d); DL <- numeric(pc)
  for (j in seq(1, pc, 2)) DL[j] <- quantile(d[, j], 0.06)
  for (j in seq_len(pc)) d[d[, j] < DL[j], j] <- 0
  expect_warning(imputeBDLs(d, dl = DL, eps = 1, method = "lm", variation = FALSE),
                 "imputeBDL")   # locale-independent (message names the replacement)
})
