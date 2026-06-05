## Deprecated aliases for the rounded-zero / below-detection-limit imputer.
## Both forward to imputeBDL() (the canonical name as of 2.6.0), preserving their
## historical default arguments so existing user code keeps working. They follow
## the package convention (see isomLR -> pivotCoord): .Deprecated() + shared
## @rdname with the replacement function.

#' @rdname imputeBDL
#' @export
imputeBDLs <- function(x, maxit = 10, eps = 0.1, method = "subPLS",
                       dl = rep(0.05, ncol(x)), variation = TRUE, nPred = NULL,
                       nComp = "boot", bruteforce = FALSE,
                       noisemethod = "residuals", noise = FALSE, R = 10,
                       correction = "normal", verbose = FALSE) {
  .Deprecated("imputeBDL")
  imputeBDL(x, maxit = maxit, eps = eps, method = method, dl = dl,
            variation = variation, nPred = nPred, nComp = nComp,
            bruteforce = bruteforce, noisemethod = noisemethod, noise = noise,
            R = R, correction = correction, verbose = verbose)
}

#' @rdname imputeBDL
#' @export
impRZilr <- function(x, maxit = 10, eps = 0.1, method = "pls",
                     dl = rep(0.05, ncol(x)), variation = FALSE, nComp = "boot",
                     bruteforce = FALSE, noisemethod = "residuals",
                     noise = FALSE, R = 10, correction = "normal",
                     verbose = FALSE) {
  .Deprecated("imputeBDL")
  imputeBDL(x, maxit = maxit, eps = eps, method = method, dl = dl,
            variation = variation, nComp = nComp, bruteforce = bruteforce,
            noisemethod = noisemethod, noise = noise, R = R,
            correction = correction, verbose = verbose)
}
