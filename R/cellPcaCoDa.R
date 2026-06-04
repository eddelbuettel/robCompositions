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
#' By default the alternating scheme is initialised from a robust MCD fit of the
#' ilr coordinates (\code{init = "mcd"}, the default in Templ 2026, underpinning
#' its breakdown results); a fast classical SVD start is available via
#' \code{init = "classical"}.  Robustness is controlled by a bounded loss
#' \code{rho}, which defaults to the Tukey biweight (95\% Gaussian efficiency),
#' as in the paper.
#'
#' The detection step uses a \emph{projection-based test} derived from the
#' propagation theorem: if raw part \eqn{j} alone is contaminated, the
#' perturbation in ilr space lies along the direction
#' \eqn{d_j = V^\top (e_j - 1/D)} (normalised to unit length).  MAD-standardised
#' ilr residuals are projected onto each \eqn{d_j}, and the projection score is
#' standardised by its correlation-adjusted scale
#' \eqn{\sigma_{p,j} = \sqrt{d_j^\top R\, d_j}}, where \eqn{R} is the correlation
#' matrix of the standardised residuals, so the standardised score is
#' approximately \eqn{N(0,1)} under the null.  A cell is flagged when this score
#' exceeds the Bonferroni-adjusted threshold \eqn{z_{1 - \alpha/(2D)}}.
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
#' \strong{Implementation choices beyond Templ (2026).}  Two conveniences are
#' not part of the published algorithm and exist only for ease of use: (i) when
#' \code{k} is \code{NULL} it is chosen by the Kaiser criterion (retain
#' eigenvalues above the average) --- the paper leaves \code{k} to the analyst,
#' so set it explicitly (e.g. from a scree plot or a clear eigenvalue gap) for
#' publication-quality results; and (ii) the ridge parameter \code{reg}
#' (default \code{0}, i.e. the paper's algorithm) adds a small diagonal term for
#' ill-conditioned high-dimensional compositions, a regime the paper lists as
#' future work.  Convergence is declared when the objective's relative change
#' falls below \code{tol}, or when a period-2 limit cycle is detected (a few
#' borderline cells oscillating under the hard-threshold reweighting); in the
#' latter case it stops at the stabilised cycle (the solution is stable up to
#' those few cells) and \code{cycled} is \code{TRUE}.
#'
#' @param x compositional data with at least three parts (a data.frame or
#'   matrix of strictly positive values whose rows are compositions, D >= 3)
#' @param k number of principal components to retain.  Cellwise detection needs
#'   a non-empty residual space, so use \code{k < D - 1}; with \code{k = D - 1}
#'   the reconstruction is exact and no cells are flagged.  If \code{NULL}
#'   (default), \code{k} is chosen by the Kaiser criterion (eigenvalues above
#'   the average) --- an implementation convenience, not part of Templ (2026);
#'   set \code{k} explicitly (e.g. from a scree plot) for publication results.
#' @param method contamination detection method: currently only
#'   \code{"adaptive"} is implemented.  A \code{"fixed"} method is planned
#'   for a future release.
#' @param rho loss function: \code{"tukey"} (default, Tukey biweight, as in
#'   Templ 2026) or \code{"huber"}.
#' @param alpha tuning constant for the loss.  If \code{NULL}, 4.685 (Tukey) or
#'   1.345 (Huber) is used, each giving 95\% Gaussian efficiency.
#' @param maxiter maximum number of alternating iterations (default: 100)
#' @param tol convergence tolerance on the relative change in the objective
#'   (default: 1e-6)
#' @param reg ridge parameter added to the diagonal of the weighted covariance
#'   for ill-conditioned high-dimensional compositions (default: 0, i.e. the
#'   published algorithm).  An implementation extension; see Details.
#' @param w0 floor weight assigned to contaminated cells (default: 0.1).
#'   Must be in (0, 1).
#' @param alpha_detect Bonferroni significance level for the cellwise
#'   contamination detection test (default: 0.01).
#' @param trace logical; if \code{TRUE}, iteration progress is printed
#' @param init initialisation of the alternating scheme: \code{"mcd"} (default)
#'   uses a robust MCD fit of the ilr coordinates (the default in Templ 2026,
#'   underpinning its breakdown results); \code{"classical"} uses a fast
#'   classical SVD start.  Falls back to \code{"classical"} (with a warning) if
#'   \code{\link[robustbase]{covMcd}} fails.
#'
#' @return An object of class \code{"cellPcaCoDa"} with components:
#'   \item{scores}{n x k score matrix in ilr space}
#'   \item{loadings}{D x k loadings in clr space}
#'   \item{eigenvalues}{the first k eigenvalues of the weighted covariance}
#'   \item{eigenvalues_all}{all eigenvalues of the final weighted covariance
#'     (length D-1), for explained-variance reporting}
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
#'   \item{init}{the initialisation actually used (\code{"mcd"} or
#'     \code{"classical"})}
#'   \item{cycled}{logical: \code{TRUE} if the algorithm stopped at a period-2
#'     limit cycle (the loadings are stable up to a few oscillating borderline
#'     cells) rather than at a fixed point; \code{converged} is then also
#'     \code{TRUE}}
#'
#' @author Matthias Templ
#' @references
#' Templ, M. (2026). cellPcaCoDa: cellwise-robust principal components for
#'   compositional data. \emph{Advances in Data Analysis and Classification}.
#'   In press.
#'
#' Templ, M. (2026). Log-ratio propagation on the simplex: a theory of cellwise
#'   contamination for compositional data. \emph{arXiv} 2605.31345.
#'
#' @note
#' \code{cellPcaCoDa} targets \emph{cellwise} contamination (individual
#' corrupted parts); use \code{\link{pcaCoDa}} when whole observations are
#' outlying (rowwise/casewise contamination).  Inspect
#' \code{colMeans(res$cellflags)} for the per-part contamination rate and
#' \code{rowSums(res$cellflags)} for the per-observation count.  The cell-weight
#' floor \code{w0} is stable on \eqn{[0.05, 0.2]}; for heavier contamination
#' (rate \eqn{\ge 0.15}) use \eqn{[0.05, 0.08]}.  The detection level
#' \code{alpha_detect} (default \code{0.01}) is Bonferroni-adjusted across the
#' \eqn{D} parts.  The input must be \strong{strictly positive} (no zeros);
#' replace rounded or structural zeros first, e.g. with \code{\link{imputeBDLs}}.
#' At least three parts are required (\eqn{D \ge 3}), and cellwise detection
#' needs \eqn{k < D-1} (at \eqn{k = D-1} the residual space is empty and no
#' cells are flagged).  The method and its defaults assume a fixed number of
#' parts with growing sample size.
#'
#' @seealso \code{\link{pcaCoDa}} for rowwise-robust PCA,
#'   \code{\link{outCoDa}} for rowwise outlier detection,
#'   \code{\link{contaminate_simplex}} to generate cellwise-contaminated data,
#'   \code{\link{pivotCoord}}, \code{\link{orthbasis}}
#' @keywords multivariate robust
#' @importFrom stats cov mad qnorm
#' @export
#' @examples
#' ## --- a 5-part composition: detection enabled (k = 2 < D-1) -------------
#' data(expenditures)
#' set.seed(123)
#' x <- expenditures
#' x[1, 2] <- x[1, 2] * 10      # contaminate one cell
#' x <- constSum(x)             # re-close to the simplex
#' res <- cellPcaCoDa(x, k = 2)
#' res                          # print method
#' head(res$cellflags)          # logical n x D contamination flags
#'
#' \donttest{
#' plot(res)                    # contamination heatmap
#' plot(res, which = 2)         # clr biplot
#'
#' ## --- paper application: GEMAS geochemistry (Templ 2026, Sec. 5) --------
#' data(gemas)
#' elements <- intersect(c("Al", "Ca", "Fe", "K", "Mg", "Mn", "Na"),
#'                       colnames(gemas))
#' xg <- gemas[, elements]
#' xg <- xg[complete.cases(xg), ]
#' xg <- constSum(xg, const = 100)               # close to 100%
#' resg <- cellPcaCoDa(xg, k = 3)                 # 7 parts, k = 3 < D-1
#' ## per-element cellwise contamination rate
#' round(sort(colMeans(resg$cellflags), decreasing = TRUE), 3)
#' plot(resg)                   # cellwise flags / clr biplot
#' }
cellPcaCoDa <- function(x, k = NULL, method = "adaptive",
                        rho = "tukey", alpha = NULL,
                        maxiter = 100, tol = 1e-6,
                        reg = 0, w0 = 0.1, alpha_detect = 0.01,
                        trace = FALSE, init = c("mcd", "classical")) {

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
  if (D < 3) {
    stop("'x' must have at least 3 parts (D >= 3); cellwise PCA is not ",
         "defined for two-part compositions")
  }
  if (n < D) {
    warning("fewer observations than parts --- results may be unstable")
  }
  method <- match.arg(method, choices = c("adaptive"))
  rho <- match.arg(rho, choices = c("huber", "tukey"))
  init <- match.arg(init)

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
  if (k >= p) {
    warning("k = D - 1: the residual space is empty, so cellwise detection ",
            "is disabled (no cells will be flagged); use k < D - 1 to enable it")
  }

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

  ## --- Step 1: Initialisation ---
  ## Default: robust MCD fit of the ilr coordinates (Templ 2026, the published
  ## default; its breakdown results rely on a robust start).  A fast classical
  ## SVD start is available via init = "classical".  MCD can fail for too few
  ## observations relative to the dimension, in which case we fall back.
  if (init == "mcd") {
    mc <- tryCatch(robustbase::covMcd(z), error = function(e) NULL)
    if (is.null(mc)) {
      warning("covMcd() failed (e.g. too few observations for the dimension); ",
              "falling back to classical initialisation")
      init <- "classical"
    }
  }
  if (init == "mcd") {
    mu <- mc$center
    zc <- sweep(z, 2, mu)
    e0 <- eigen(mc$cov, symmetric = TRUE)
    P  <- e0$vectors[, seq_len(k), drop = FALSE]   # p x k loadings in ilr space
    scores <- zc %*% P                             # n x k scores
  } else {
    mu <- colMeans(z)
    zc <- sweep(z, 2, mu)
    sv <- svd(zc, nu = k, nv = k)
    P  <- sv$v                        # p x k loadings in ilr space
    scores <- zc %*% P                # n x k scores
  }

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
  obj_old  <- Inf
  obj_old2 <- Inf
  converged <- FALSE
  cycled    <- FALSE

  for (iter in seq_len(maxiter)) {

    ## (a) Compute residuals in ilr space
    fitted <- tcrossprod(scores, P)             # n x p
    resid  <- zc - fitted
    sigma  <- apply(resid, 2, mad)
    sigma[sigma < 1e-10] <- 1e-10
    resid_std <- sweep(resid, 2, sigma, "/")

    ## (b) Detect contaminated cells in RAW PART space
    ## Projection-based detection (Templ 2026, Eq. 17): project the MAD-
    ## standardised residuals onto each unit direction d_j, then standardise the
    ## projection by its correlation-adjusted scale
    ## sigma_{p,j} = sqrt(d_j^T R d_j), where R is the correlation matrix of the
    ## standardised residuals, so the score is approximately N(0,1) under the
    ## null.  A cell is flagged against a Bonferroni-adjusted normal quantile.
    proj_scores <- resid_std %*% proj_dirs        # n x D ; p_ij = d_j^T r_i
    R_corr <- stats::cor(resid_std)               # p x p correlation matrix
    R_corr[!is.finite(R_corr)] <- 0               # guard zero-variance/NA cols
    diag(R_corr) <- 1
    sigma_p <- sqrt(pmax(colSums((R_corr %*% proj_dirs) * proj_dirs),
                         .Machine$double.eps))    # length D: sqrt(d_j^T R d_j)
    proj_std <- sweep(proj_scores, 2, sigma_p, "/")  # standardised scores
    thresh <- qnorm(1 - alpha_detect / (2 * D))   # Bonferroni over D parts
    cellflags <- abs(proj_std) > thresh

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
    col_w[col_w < 1e-10] <- 1e-10                # guard against zero weights
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
    if (is.nan(obj_new) || is.na(obj_new)) obj_new <- Inf
    if (trace) {
      message(sprintf("Iteration %d: objective = %.6f, flagged = %d / %d cells",
                      iter, obj_new, sum(cellflags), n * D))
    }
    ## (i) standard convergence: consecutive objective change below tolerance
    if (is.finite(obj_old) && is.finite(obj_new) &&
        abs(obj_old - obj_new) / max(1, abs(obj_old)) < tol) {
      converged <- TRUE
      break
    }
    ## (ii) period-2 limit cycle: the objective matches the value two iterations
    ## back while the consecutive change does not.  The hard-threshold
    ## reweighting is oscillating over a few borderline cells; treat as
    ## stabilised and return the better cycled iterate.
    if (is.finite(obj_old2) && is.finite(obj_new) &&
        abs(obj_old2 - obj_new) / max(1, abs(obj_old2)) < tol) {
      converged <- TRUE
      cycled    <- TRUE
      break
    }
    obj_old2 <- obj_old
    obj_old  <- obj_new
  }

  ## On a detected limit cycle we keep the current (stabilised) iterate: the
  ## objective is not monotone in solution quality here (it is smallest for the
  ## over-flagged first iterate), so a min-objective rule would be misleading.

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
    rho          = rho,
    init         = init,
    cycled       = cycled
  )
  class(res) <- "cellPcaCoDa"
  invisible(res)
}
