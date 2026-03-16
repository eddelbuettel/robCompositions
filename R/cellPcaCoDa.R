#' Cellwise-robust PCA for compositional data
#'
#' Performs principal component analysis for compositional data that is robust
#' against cellwise contamination in raw parts.  The method accounts for the
#' dependent contamination structure induced by closure on the simplex.
#'
#' Compositional data are first expressed in isometric logratio (pivot)
#' coordinates via \code{\link{pivotCoord}}.  An alternating optimisation
#' scheme then jointly estimates cell-level contamination flags in the original
#' \emph{raw-part} space, maps them to coordinate-level weights through the
#' contrast matrix, and computes a weighted PCA in the ilr space.  Loadings
#' and the robust centre are back-transformed to the clr space for
#' interpretation, following the same convention as \code{\link{pcaCoDa}}.
#'
#' The detection step uses a \emph{projection-based test} derived from the
#' propagation theorem: if raw part \eqn{j} alone is contaminated, the
#' perturbation in ilr space lies along the direction
#' \eqn{d_j = V^\top (e_j - 1/D)}.  Standardised residuals are projected
#' onto each \eqn{d_j} and compared with a Bonferroni-adjusted threshold.
#'
#' The weight mapping from raw-part flags to ilr weights is fully vectorised
#' via a matrix multiply (no triple-nested loop).  For each observation
#' \eqn{i} and ilr coordinate \eqn{l} the weight is
#' \deqn{w_{il} = \max\Bigl(w_0,\;
#'   \exp\bigl(\sum_j c_{ij}\,\log\bigl(1-(1-w_0)\,|V_{jl}|/\max_j |V_{jl}|\bigr)\bigr)\Bigr)}
#' where \eqn{c_{ij}} are the binary contamination flags (1 = flagged).
#'
#' The weighted covariance uses \emph{per-coordinate} normalisation:
#' entry \eqn{(l, m)} of the covariance is normalised by
#' \eqn{\sqrt{\sum_i w_{il} \cdot \sum_i w_{im}}} so that coordinates with
#' different total weights are treated correctly.
#'
#' The default number of components \code{k} is chosen by the Kaiser
#' criterion (retain eigenvalues above the average).  This is a simple
#' heuristic; alternatives include parallel analysis (Horn, 1965) and
#' BIC-based selection.  Users should set \code{k} explicitly for
#' publication-quality results.
#'
#' @param x compositional data (data.frame or matrix with strictly positive
#'   values whose rows represent compositions)
#' @param k number of principal components to retain.  If \code{NULL}
#'   (default), the Kaiser criterion is used (eigenvalues above the average).
#' @param method contamination detection method: currently only
#'   \code{"adaptive"} is implemented.  A \code{"fixed"} method is planned
#'   for a future release.
#' @param rho loss function: \code{"huber"} (default) or \code{"tukey"}
#' @param alpha tuning parameter for the loss function.  If \code{NULL},
#'   1.345 (Huber, 95\% efficiency) or 4.685 (Tukey) is used.
#' @param maxiter maximum number of alternating iterations (default: 100)
#' @param tol convergence tolerance on relative change in the objective
#'   (default: 1e-6)
#' @param reg regularisation parameter added to the diagonal of the
#'   weighted covariance for high-dimensional compositions (default: 0)
#' @param w0 floor weight assigned to contaminated cells (default: 0.1).
#'   Must be in (0, 1).
#' @param alpha_detect Bonferroni significance level for the cellwise
#'   contamination detection test (default: 0.01).
#' @param trace logical; if \code{TRUE}, iteration progress is printed
#'
#' @return An object of class \code{"cellPcaCoDa"} with components:
#'   \item{scores}{n x k score matrix in ilr space}
#'   \item{loadings}{D x k loadings in clr space}
#'   \item{eigenvalues}{the first k eigenvalues of the weighted covariance}
#'   \item{cellflags}{n x D logical matrix of contamination flags in
#'     raw-part space}
#'   \item{cellweights}{n x (D-1) matrix of cell weights in ilr space}
#'   \item{center}{robust centre in ilr space (length D-1)}
#'   \item{center_clr}{robust centre in clr space (length D)}
#'   \item{converged}{logical: did the algorithm converge within
#'     \code{maxiter} iterations?}
#'   \item{iterations}{number of iterations actually used}
#'   \item{method}{the detection method used}
#'   \item{k}{number of retained components}
#'   \item{rho}{the loss function name}
#'
#' @author Matthias Templ
#' @references Templ, M. (2026).
#'   Cellwise contamination on the simplex: theory and robust PCA for
#'   compositions.
#'
#' @seealso \code{\link{pcaCoDa}} for rowwise-robust PCA,
#'   \code{\link{outCoDa}} for rowwise outlier detection,
#'   \code{\link{pivotCoord}}, \code{\link{orthbasis}}
#' @keywords multivariate robust
#' @importFrom stats cov mad qnorm
#' @export
#' @examples
#' data(arcticLake)
#'
#' ## introduce cellwise contamination
#' set.seed(123)
#' x <- arcticLake
#' x[1, 2] <- x[1, 2] * 10   # contaminate one cell
#' x <- constSum(x)            # re-close
#' res <- cellPcaCoDa(x, k = 2)
#' res                          # print method
#' res$cellflags                # should flag row 1, part 2
#' plot(res)                    # contamination heatmap
#' plot(res, which = 2)         # biplot
#'
#' ## clean data, classical default
#' data(expenditures)
#' res2 <- cellPcaCoDa(expenditures)
#' res2
cellPcaCoDa <- function(x, k = NULL, method = "adaptive",
                        rho = "huber", alpha = NULL,
                        maxiter = 100, tol = 1e-6,
                        reg = 0, w0 = 0.1, alpha_detect = 0.01,
                        trace = FALSE) {

  ## --- Input validation ---
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data.frame")
  }
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("'x' must be numeric")
  }
  n <- nrow(x)
  D <- ncol(x)
  if (any(x <= 0, na.rm = TRUE)) {
    stop("all values in 'x' must be strictly positive")
  }
  if (anyNA(x)) {
    stop("missing values in 'x' are not supported; consider imputation first")
  }
  if (n < D) {
    warning("fewer observations than parts --- results may be unstable")
  }
  method <- match.arg(method, choices = c("adaptive"))
  rho <- match.arg(rho, choices = c("huber", "tukey"))

  ## --- Step 0: Transform to ilr coordinates ---
  z <- as.matrix(pivotCoord(x))   # n x (D-1)
  p <- D - 1

  ## Contrast matrix V: D x (D-1), maps clr -> ilr
  V <- orthbasis(D)$V             # D x (D-1)

  ## Default k: Kaiser criterion (retain components with above-average
  ## eigenvalues of the classical covariance).  Alternatives include
  ## parallel analysis (Horn 1965) and BIC-based selection.
  ## Users should set k explicitly for publication results.
  if (is.null(k)) {
    ev <- eigen(cov(z))$values
    k <- max(1L, sum(ev > mean(ev)))
  }
  k <- min(k, p)

  ## Default tuning constant (chosen for 95% asymptotic efficiency at the
  ## normal model)
  if (is.null(alpha)) {
    alpha <- if (rho == "huber") 1.345 else 4.685
  }

  ## --- Loss and psi functions (vectorised) ---
  rho_fun <- switch(rho,
    huber = function(u) ifelse(abs(u) <= alpha,
                               u^2 / 2,
                               alpha * (abs(u) - alpha / 2)),
    tukey = function(u) ifelse(abs(u) <= alpha,
                               (alpha^2 / 6) * (1 - (1 - (u / alpha)^2)^3),
                               alpha^2 / 6),
    stop("'rho' must be \"huber\" or \"tukey\"")
  )
  psi_fun <- switch(rho,
    huber = function(u) ifelse(abs(u) <= alpha,
                               u,
                               alpha * sign(u)),
    tukey = function(u) ifelse(abs(u) <= alpha,
                               u * (1 - (u / alpha)^2)^2,
                               0)
  )

  ## --- Step 1: Classical PCA initialisation ---
  mu <- colMeans(z)
  zc <- sweep(z, 2, mu)
  sv <- svd(zc, nu = k, nv = k)
  P  <- sv$v                        # p x k loadings in ilr space
  scores <- zc %*% P                # n x k scores

  ## Initialise cell flags: no contamination
  cellflags <- matrix(FALSE, nrow = n, ncol = D)
  W <- matrix(1, nrow = n, ncol = p)  # ilr cell weights

  ## --- Precompute projection directions for detection step ---
  ## d_j = V^T (e_j - 1/D), normalised to unit length.
  ## These are the theoretical perturbation directions from the
  ## propagation theorem (one per raw part j = 1, ..., D).
  proj_dirs <- matrix(0, nrow = p, ncol = D)
  for (j in seq_len(D)) {
    d_j <- as.numeric(crossprod(V, diag(D)[, j] - 1 / D))
    proj_dirs[, j] <- d_j / sqrt(sum(d_j^2))
  }

  ## Precompute max-normalised |V| for vectorised weight mapping
  V_abs <- abs(V)                             # D x p
  V_colmax <- apply(V_abs, 2, max)
  V_colmax[V_colmax < .Machine$double.eps] <- 1
  V_norm <- sweep(V_abs, 2, V_colmax, "/")   # D x p, columns max-normalised

  ## --- Step 2: Alternating optimisation ---
  obj_old <- Inf
  converged <- FALSE

  for (iter in seq_len(maxiter)) {

    ## (a) Compute residuals in ilr space
    fitted <- tcrossprod(scores, P)             # n x p
    resid  <- zc - fitted
    sigma  <- apply(resid, 2, mad)
    sigma[sigma < 1e-10] <- 1e-10
    resid_std <- sweep(resid, 2, sigma, "/")

    ## (b) Detect contaminated cells in RAW PART space
    ## Projection-based detection: project standardised residuals onto
    ## each d_j and apply Bonferroni-corrected threshold.
    proj_scores <- resid_std %*% proj_dirs     # n x D
    thresh <- qnorm(1 - alpha_detect / (2 * D)) # Bonferroni-adjusted
    cellflags <- abs(proj_scores) > thresh

    ## (c) Map raw-part flags to ilr weights --- vectorised (Issue 32 fix)
    ## Uses matrix multiply instead of triple-nested loop: O(n*D*(D-1)) as
    ## a single matmul rather than element-wise branching.
    flag_matrix    <- cellflags * 1.0                      # n x D, numeric 0/1
    log_reduction  <- flag_matrix %*% log(1 - (1 - w0) * V_norm)  # n x p
    W <- pmax(exp(log_reduction), w0)

    ## (d) Weighted PCA update
    ## Weighted centering
    mu <- colSums(W * z) / colSums(W)
    zc <- sweep(z, 2, mu)

    ## Weighted covariance --- per-coordinate normalisation (Issue 6 fix)
    ## Each (l, m) entry of S is normalised by sqrt(sum(W[,l]) * sum(W[,m]))
    ## so that coordinates with different total weights are treated correctly.
    Wz <- sqrt(W) * zc
    col_w   <- colSums(W)                        # length-p vector
    normmat <- sqrt(outer(col_w, col_w))          # p x p normalisation matrix
    S <- crossprod(Wz) / normmat
    if (reg > 0) {
      S <- S + reg * diag(p)
    }

    ## Eigen-decomposition
    eig <- eigen(S, symmetric = TRUE)
    if (any(eig$values < 0)) {
      warning("weighted covariance has negative eigenvalues; consider increasing 'reg'")
    }
    P   <- eig$vectors[, seq_len(k), drop = FALSE]
    eigenvalues <- eig$values[seq_len(k)]
    scores <- zc %*% P

    ## (e) Check convergence on weighted objective
    ## Recompute residuals with updated PCA for consistent objective
    fitted_new    <- tcrossprod(scores, P)
    resid_new     <- zc - fitted_new
    sigma_new     <- apply(resid_new, 2, mad)
    sigma_new[sigma_new < 1e-10] <- 1e-10
    resid_std_new <- sweep(resid_new, 2, sigma_new, "/")
    obj_new <- sum(W * rho_fun(resid_std_new))
    if (trace) {
      message(sprintf("Iteration %d: objective = %.6f, flagged = %d / %d cells",
                      iter, obj_new, sum(cellflags), n * D))
    }
    if (abs(obj_old - obj_new) / max(1, abs(obj_old)) < tol) {
      converged <- TRUE
      break
    }
    obj_old <- obj_new
  }

  ## --- Step 3: Back-transform to clr space ---
  loadings_clr <- V %*% P                     # D x k
  center_clr   <- as.numeric(V %*% mu)        # length D

  if (!is.null(colnames(x))) {
    rownames(loadings_clr) <- colnames(x)
    colnames(cellflags)    <- colnames(x)
    names(center_clr)      <- colnames(x)
  }
  colnames(loadings_clr) <- paste0("Comp.", seq_len(k))

  ## Store all eigenvalues for explained-variance reporting
  eigenvalues_all <- eig$values

  ## --- Return ---
  res <- list(
    scores       = scores,
    loadings     = loadings_clr,
    eigenvalues  = eigenvalues,
    eigenvalues_all = eigenvalues_all,
    cellflags    = cellflags,
    cellweights  = W,
    center       = mu,
    center_clr   = center_clr,
    converged    = converged,
    iterations   = iter,
    method       = method,
    k            = k,
    rho          = rho
  )
  class(res) <- "cellPcaCoDa"
  invisible(res)
}
