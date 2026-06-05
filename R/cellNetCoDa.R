## cellNetCoDa: cellwise-robust network (Gaussian graphical model) estimation
## for compositional data. Public entry point cellNetCoDa() plus the internal
## pipeline (zero treatment -> cellwise detection -> cell-cleaned covariance ->
## StARS penalty selection -> MB / graphical-lasso back end -> graph). The S3
## print/summary/plot methods live in R/plot.cellNetCoDa.R.

#' Assemble the taxa-level graph object from a selected adjacency.
#' @noRd
to_taxa_graph <- function(adj, taxa = NULL, weights = NULL, extra = list()) {
  D <- nrow(adj)
  if (is.null(taxa)) taxa <- paste0("T", seq_len(D))
  dimnames(adj) <- list(taxa, taxa)
  out <- c(list(adj = adj, taxa = taxa, n_edges = sum(adj[upper.tri(adj)]),
                weights = weights), extra)
  class(out) <- "cellNetCoDa"
  out
}

#' Cellwise-robust network estimation for compositional data
#'
#' Estimates a sparse conditional-independence (Gaussian graphical) network
#' among the parts of a composition that is robust against \emph{cellwise}
#' contamination --- individual corrupted parts --- on top of the closure
#' constraint of the simplex.
#'
#' The estimator chains four steps.  (i) Count zeros are replaced by a
#' Bayesian-multiplicative replacement and the rows are closed to the simplex.
#' (ii) Cellwise outliers are flagged by a masking-resistant conditional-pivot
#' test together with the projection flags of \code{\link{cellPcaCoDa}} and
#' turned into per-cell weights; a cell-cleaned centred-logratio (clr)
#' covariance is then formed by shrinking each flagged cell toward the robust
#' clr centre,
#' \deqn{\tilde{C}_{ij} = w_{ij}\,C_{ij} + (1 - w_{ij})\,\mu_j,}
#' where \eqn{C} is the clr matrix, \eqn{w_{ij}\in[0,1]} the cell weight and
#' \eqn{\mu} the robust clr centre.  (iii) A sparse graph is read off the
#' cleaned covariance, either by cellwise-robust Meinshausen--Buhlmann
#' neighbourhood selection (\code{method = "mb"}, the default) or by the
#' graphical lasso on the clr covariance (\code{method = "glasso"}).  (iv) The
#' sparsity penalty is chosen by StARS.
#'
#' \strong{Penalty selection (StARS).}  For each penalty on the supplied grid
#' the graph is refitted on \code{B} random subsamples (fraction
#' \code{sub_ratio}); the edge-selection frequency \eqn{\theta} gives the edge
#' instability \eqn{2\theta(1-\theta)}, averaged over edges.  Following Liu et
#' al. (2010) the instability is monotonised and the \emph{densest} graph whose
#' instability stays at or below \code{thresh} (the smallest stable penalty) is
#' returned.
#'
#' \strong{Back ends.}  \code{"mb"} regresses each part's clr value on the
#' cell-cleaned clr of the remaining parts by a weighted lasso and combines the
#' neighbourhoods by \code{rule}; \code{"glasso"} penalises the clr covariance
#' directly and reads the graph from the exact zeros of the estimated precision
#' matrix (the SPIEC-EASI convention).  Both share the same subsampling, grid
#' and instability convention, so penalties are comparable across back ends.
#'
#' @param X an \eqn{n \times D} matrix or data.frame of compositional parts
#'   (columns are parts).  Counts or proportions are accepted; count zeros are
#'   pre-treated by a Bayesian-multiplicative replacement
#'   (\code{\link[zCompositions]{cmultRepl}}, method \code{"GBM"}) and the rows
#'   closed to the simplex internally.  At least three parts (\eqn{D \ge 3}).
#' @param method back end for the sparse network estimate: \code{"mb"}
#'   (default) for cellwise-robust Meinshausen--Buhlmann neighbourhood
#'   selection in clr space, or \code{"glasso"} for the graphical lasso on the
#'   clr covariance.
#' @param rho_grid numeric vector of graphical-lasso penalties (decreasing
#'   order recommended).  \strong{Required when \code{method = "glasso"}};
#'   ignored otherwise.
#' @param lambda_grid numeric vector of neighbourhood-selection penalties
#'   (decreasing).  \strong{Required when \code{method = "mb"}}; ignored
#'   otherwise.
#' @param B number of subsamples used by the StARS penalty selection
#'   (default: 20).
#' @param sub_ratio subsample fraction for StARS (default: 0.7).
#' @param thresh StARS edge-instability threshold (default: 0.05).  The densest
#'   graph whose monotonised instability stays at or below \code{thresh} is
#'   selected.
#' @param k number of principal components passed to the cellwise detector
#'   \code{\link{cellPcaCoDa}}.  If \code{NULL} (default) it is chosen there by
#'   the Kaiser criterion; set it explicitly (\code{k < D - 1}) for
#'   reproducible results.
#' @param rule neighbourhood-combination rule for \code{method = "mb"}:
#'   \code{"and"} (default; an edge requires both nodes to select each other)
#'   or \code{"or"}.  Ignored for \code{"glasso"}.
#' @param taxa optional character vector of node (part) names for the returned
#'   graph; defaults to \code{colnames(X)}, or \code{T1, \dots, TD} when
#'   \code{X} is unnamed.
#' @param \dots further arguments for the cellwise detector, e.g.
#'   \code{grading} (\code{"marginal"}, \code{"projection"} or
#'   \code{"isotropic"}), the flag floor \code{w0}, the detection cutoff
#'   \code{k_thr}, or the number of de-masking passes \code{iter}; remaining
#'   arguments are passed on to \code{\link{cellPcaCoDa}}.  These are used both
#'   in the StARS subsampling and in the final fit.
#'
#' @return An object of class \code{"cellNetCoDa"} with components:
#'   \item{adj}{\eqn{D \times D} symmetric 0/1 adjacency matrix with zero
#'     diagonal and \code{dimnames} given by \code{taxa}: the estimated
#'     conditional-independence graph of the parts}
#'   \item{taxa}{character vector of node (part) names}
#'   \item{n_edges}{number of edges (the upper-triangle sum of \code{adj})}
#'   \item{weights}{\eqn{n \times D} matrix of cellwise weights in
#'     \eqn{[0, 1]} (1 = clean, smaller = more strongly down-weighted)}
#'   \item{penalty}{the StARS-selected penalty (\code{rho} for glasso,
#'     \code{lambda} for MB)}
#'   \item{instability}{the StARS instability curve over the penalty grid}
#'   \item{grid}{the penalty grid the StARS curve was evaluated on}
#'   \item{thresh}{the StARS instability threshold used for selection}
#'   \item{method}{the back end used, \code{"mb"} or \code{"glasso"}}
#'   \item{flags}{\eqn{n \times D} logical matrix of cellwise contamination
#'     flags}
#'
#' @author Matthias Templ
#' @references
#' Liu, H., Roeder, K. and Wasserman, L. (2010). Stability approach to
#'   regularization selection (StARS) for high-dimensional graphical models.
#'   \emph{Advances in Neural Information Processing Systems} 23.
#'
#' Meinshausen, N. and Buhlmann, P. (2006). High-dimensional graphs and
#'   variable selection with the lasso. \emph{The Annals of Statistics}
#'   34(3), 1436--1462.
#'
#' Friedman, J., Hastie, T. and Tibshirani, R. (2008). Sparse inverse
#'   covariance estimation with the graphical lasso. \emph{Biostatistics}
#'   9(3), 432--441.
#'
#' Kurtz, Z. D., Muller, C. L., Miraldi, E. R., Littman, D. R., Blaser, M. J.
#'   and Bonneau, R. A. (2015). Sparse and compositionally robust inference of
#'   microbial ecological networks. \emph{PLoS Computational Biology}
#'   11(5), e1004226.
#'
#' Raymaekers, J. and Rousseeuw, P. J. (2024). The cellwise minimum covariance
#'   determinant estimator. \emph{Journal of the American Statistical
#'   Association} 119(548), 2610--2621.
#'
#' Martin-Fernandez, J. A., Hron, K., Templ, M., Filzmoser, P. and
#'   Palarea-Albaladejo, J. (2015). Bayesian-multiplicative treatment of count
#'   zeros in compositional data sets. \emph{Statistical Modelling}
#'   15(2), 134--158.
#'
#' @note
#' The input must be compositional (positive parts); count zeros are handled
#' internally, but structural / below-detection zeros are better addressed
#' first (e.g. \code{\link{imputeBDL}}).  At least three parts are required.
#' \code{method = "mb"} needs \code{lambda_grid} and \code{method = "glasso"}
#' needs \code{rho_grid}; supplying neither is an error.  Although the nodes
#' are labelled \code{taxa} (the method originates in microbiome network
#' analysis), the parts may be any compositional variables.  StARS is
#' stochastic --- use \code{set.seed} for reproducible penalty selection.
#'
#' @seealso \code{\link{cellPcaCoDa}} for the cellwise detector it builds on,
#'   \code{\link{pcaCoDa}} and \code{\link{outCoDa}} for rowwise robustness,
#'   \code{\link{pivotCoord}}, \code{\link{orthbasis}}
#' @keywords multivariate robust
#' @export
#' @examples
#' \donttest{
#' data(expenditures)            # 5-part compositional budget data
#' set.seed(1)
#' ## cellwise-robust network via MB neighbourhood selection (needs 'glmnet')
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   net <- cellNetCoDa(expenditures, method = "mb",
#'                      lambda_grid = c(0.4, 0.3, 0.2, 0.1, 0.05),
#'                      B = 10, k = 2)
#'   net                         # print method: back end, edges, penalty
#'   net$adj                     # estimated 0/1 adjacency among the parts
#' }
#'
#' ## graphical-lasso back end, SPIEC-EASI style (needs 'glasso')
#' if (requireNamespace("glasso", quietly = TRUE)) {
#'   net_gl <- cellNetCoDa(expenditures, method = "glasso",
#'                         rho_grid = c(0.4, 0.3, 0.2, 0.1, 0.05),
#'                         B = 10, k = 2)
#'   net_gl$n_edges
#' }
#' }
cellNetCoDa <- function(X, method = c("mb", "glasso"),
                        rho_grid = NULL, lambda_grid = NULL,
                        B = 20, sub_ratio = 0.7, thresh = 0.05, k = NULL,
                        rule = "and", taxa = colnames(X), ...) {
  method <- match.arg(method)
  ## glmnet (method "mb") and glasso (method "glasso") are Suggests; check the
  ## one needed here so the user gets a clear message before any subsampling.
  pkg <- if (method == "glasso") "glasso" else "glmnet"
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("cellNetCoDa(method = \"", method, "\") requires the '", pkg,
         "' package; install it with install.packages(\"", pkg, "\").",
         call. = FALSE)
  }
  X0 <- treat_zeros(X)
  sel <- stars_select(X0, method = method, rho_grid = rho_grid, lambda_grid = lambda_grid,
                      B = B, sub_ratio = sub_ratio, thresh = thresh, k = k, rule = rule, ...)
  detF <- detect_cells(X0, k = k, ...)
  to_taxa_graph(sel$adj, taxa = taxa, weights = detF$weights,
                extra = list(penalty = sel$penalty, instability = sel$instability,
                             grid = sel$grid, thresh = thresh,
                             method = method, flags = detF$flags))
}


## ---- helpers.R ----

# Coordinate helpers for cellNetCoDa. clr is the single source of basis truth;
# ilr = clr %*% V, clr = ilr %*% t(V), with V = orthbasis(D)$V (D x (D-1), V'V = I).

clr_mat <- function(X) {
  X <- as.matrix(X)
  L <- log(X)
  L - rowMeans(L)                       # n x D, rows sum to 0
}

clr_inv_mat <- function(C) {
  E <- exp(as.matrix(C)); E / rowSums(E)
}

get_V <- function(D) orthbasis(D)$V   # D x (D-1), orthonormal

ilr_from_clr <- function(C, V) as.matrix(C) %*% V       # n x (D-1)
clr_from_ilr <- function(Z, V) as.matrix(Z) %*% t(V)    # n x D (centred)

close_comp <- function(X) constSum(as.matrix(X))   # closure to sum 1

# Tukey bisquare weight in [0,1] for standardised residual u; w=0 beyond c.
tukey_w <- function(u, c = 3) ifelse(abs(u) < c, (1 - (u / c)^2)^2, 0)

## ---- treat_zeros.R ----

#' Replace count zeros by a Bayesian-multiplicative replacement, then close.
#'
#' @param X numeric matrix or data.frame of counts or proportions (parts in
#'   columns).
#' @param method replacement method passed to
#'   \code{zCompositions::cmultRepl()} (default \code{"GBM"}).
#' @param \dots further arguments for \code{zCompositions::cmultRepl()}.
#' @return the zero-replaced data closed to the simplex.
#' @noRd
treat_zeros <- function(X, method = "GBM", ...) {
  X <- as.matrix(X)
  if (any(X == 0)) {
    Xt <- zCompositions::cmultRepl(X, label = 0, method = method,
                                   output = "prop", suppress.print = TRUE, ...)
    return(close_comp(as.matrix(Xt)))
  }
  close_comp(X)
}

## ---- robust_ilr_cov.R ----

#' Cell-cleaned ilr covariance from detector weights.
#'
#' Shrinks each clr cell toward the robust clr centre by its detector weight
#' (\eqn{w C + (1-w)\mu}), maps the cleaned data to ilr coordinates, and returns
#' the covariance of the cleaned coordinates together with the raw and cleaned
#' coordinate matrices.
#'
#' @param X compositional data (closed internally).
#' @param det output of \code{detect_cells()}; supplies the cell \code{weights}
#'   and the robust clr centre \code{mu_clr}.
#' @param V optional clr-to-ilr contrast matrix; defaults to
#'   \code{orthbasis(D)$V}.
#' @return list with \code{Sigma} (ilr covariance of the cleaned data),
#'   \code{Z} and \code{Zc} (raw and cleaned ilr coordinates), \code{Cc} (the
#'   cleaned clr matrix) and \code{V}.
#' @noRd
robust_ilr_cov <- function(X, det, V = NULL) {
  X <- close_comp(as.matrix(X)); D <- ncol(X)
  if (is.null(V)) V <- get_V(D)
  W <- det$weights; mu <- det$mu_clr
  C <- clr_mat(X)
  Cc <- W * C + (1 - W) * matrix(mu, nrow(C), D, byrow = TRUE)   # cell-wise shrink to robust centre
  Z  <- C  %*% V
  Zc <- Cc %*% V
  list(Sigma = stats::cov(Zc), Z = Z, Zc = Zc, Cc = Cc, V = V)
}

## ---- precision_mb.R ----

#' Cellwise-robust Meinshausen-Buhlmann neighbourhood selection in clr space.
#'
#' For each part the clr response is regressed on the cell-cleaned clr of the
#' remaining parts by a weighted lasso (cell weights from \code{detect_cells});
#' the per-part neighbourhoods are combined into a symmetric adjacency by
#' \code{rule}.
#'
#' @param X compositional data (closed internally).
#' @param det output of \code{detect_cells()} (supplies the cell weights).
#' @param rc output of \code{robust_ilr_cov()} (supplies the cleaned clr
#'   predictors \code{rc$Cc}).
#' @param lambda optional common penalty; if \code{NULL}, chosen per part by
#'   5-fold cross-validation (\code{lambda.1se}).
#' @param rule neighbourhood-combination rule, \code{"and"} (default) or
#'   \code{"or"}.
#' @param \dots currently unused.
#' @return list with \code{adj} (D x D 0/1 adjacency, zero diagonal) and
#'   \code{Beta} (the neighbourhood coefficient matrix).
#' @noRd
precision_mb <- function(X, det, rc, lambda = NULL, rule = c("and", "or"), ...) {
  rule <- match.arg(rule)
  X <- close_comp(as.matrix(X)); D <- ncol(X)
  C  <- clr_mat(X)                      # raw clr responses
  Cc <- rc$Cc                           # cell-cleaned clr predictors
  W  <- det$weights
  Beta <- matrix(0, D, D)
  for (a in seq_len(D)) {
    y  <- C[, a]; Xp <- Cc[, -a, drop = FALSE]; wts <- W[, a]
    if (is.null(lambda)) {
      fit <- glmnet::cv.glmnet(Xp, y, weights = wts, nfolds = 5, standardize = TRUE)
      b <- as.numeric(stats::coef(fit, s = "lambda.1se"))[-1]
    } else {
      fit <- glmnet::glmnet(Xp, y, weights = wts, lambda = lambda, standardize = TRUE)
      b <- as.numeric(stats::coef(fit))[-1]
    }
    Beta[a, -a] <- b
  }
  Aedge <- (Beta != 0)
  adj <- if (rule == "and") (Aedge & t(Aedge)) else (Aedge | t(Aedge))
  diag(adj) <- FALSE
  list(adj = adj * 1L, Beta = Beta)
}

## ---- precision_glasso.R ----

#' Graphical lasso on the clr covariance -> taxa graph (the SPIEC-EASI convention).
#'
#' Runs the graphical lasso on the cell-cleaned D x D clr covariance and reads
#' the graph directly from the exact-zero sparsity pattern of the L1-penalised
#' precision matrix (off-diagonal zeros are the absent edges; no basis mapping
#' and no magnitude threshold).
#'
#' @param rc output of \code{robust_ilr_cov()} (uses \code{rc$Cc}, the cleaned
#'   clr matrix).
#' @param rho graphical-lasso penalty (on the clr correlation scale).
#' @param ridge_frac ridge fraction for the singular-covariance convergence
#'   fallback: if the fit fails, retry on
#'   \code{S_clr + ridge_frac * mean(diag(S_clr)) * I} (default 1e-3; no ridge
#'   is used when the fit converges).
#' @return list with \code{adj} (D x D 0/1 adjacency, zero diagonal),
#'   \code{Omega_clr} (the precision matrix) and \code{ridge} (the ridge
#'   actually used, 0 if none).
#' @noRd
precision_glasso <- function(rc, rho = 0.1, ridge_frac = 1e-3) {
  # Cleaned clr covariance. Use the centring projector to keep it in the clr hyperplane:
  # cov(rc$Cc) can drift off the hyperplane when the robust centre rc$mu_clr is not exactly
  # clr-centred, so we recentre. (Equivalently rc$V %*% rc$Sigma %*% t(rc$V).)
  Cc <- sweep(rc$Cc, 1, rowMeans(rc$Cc), "-")    # re-impose row-sum-zero (clr hyperplane)
  S_clr <- stats::cov(Cc)
  glasso_read_clr(S_clr, rho = rho, ridge_frac = ridge_frac)
}

#' Read a taxa graph from a clr covariance via the graphical lasso.
#'
#' Standardises the (possibly rank-deficient) clr covariance to a correlation
#' so the penalty is comparable across coordinate pairs, runs the graphical
#' lasso, and returns the exact-zero precision support as the adjacency. The
#' edge support is scale-equivariant, so the correlation standardisation does
#' not change which edges are zero at a given relative penalty. A small ridge
#' is added only as a convergence fallback for singular input.
#'
#' @param S_clr D x D clr covariance (may be singular, rank D-1).
#' @param rho graphical-lasso penalty on the correlation scale.
#' @param ridge_frac ridge fraction for the singular-matrix convergence
#'   fallback (default 1e-3).
#' @return list with \code{adj}, \code{Omega_clr} and \code{ridge}, as in
#'   \code{precision_glasso}.
#' @noRd
glasso_read_clr <- function(S_clr, rho, ridge_frac = 1e-3) {
  D <- nrow(S_clr)
  S_clr <- (S_clr + t(S_clr)) / 2                       # symmetrise numerical noise
  fit_corr <- function(S) {
    d <- sqrt(diag(S)); d[d < 1e-12] <- 1e-12
    R <- S / tcrossprod(d)                              # clr correlation
    R <- (R + t(R)) / 2; diag(R) <- 1
    gl <- glasso::glasso(R, rho = rho, penalize.diagonal = FALSE)
    list(wi = gl$wi, niter = gl$niter)
  }
  ridge_used <- 0
  res <- tryCatch(fit_corr(S_clr),
                  error = function(e) NULL, warning = function(w) NULL)
  if (is.null(res)) {
    # Singular-matrix fallback: tiny SPIEC-EASI-style ridge, then retry.
    ridge_used <- ridge_frac * mean(diag(S_clr))
    res <- fit_corr(S_clr + ridge_used * diag(D))
  }
  wi  <- (res$wi + t(res$wi)) / 2                        # symmetrise (glasso wi is not exactly symm)
  adj <- (wi != 0); diag(adj) <- FALSE
  list(adj = adj * 1L, Omega_clr = wi, ridge = ridge_used)
}

## ---- stars_select.R ----

#' Generic StARS penalty selection for an arbitrary edge estimator.
#'
#' Runs StARS on any estimator expressed as a closure
#' \code{(Xsub, penalty) -> D x D 0/1 adjacency (zero diagonal)}, so different
#' back ends can share the same subsampling, grid and instability convention.
#'
#' @param X n x D composition matrix (closed internally).
#' @param estimator function(Xs, penalty) -> D x D 0/1 adjacency (zero diagonal).
#' @param grid numeric vector of penalty values (decreasing order recommended).
#' @param B number of subsamples (default 20).
#' @param sub_ratio subsample fraction (default 0.7).
#' @param thresh instability threshold (default 0.05).
#' @return list with \code{adj}, \code{penalty}, \code{instability}, \code{grid}.
#' @noRd
stars_generic <- function(X, estimator, grid, B = 20, sub_ratio = 0.7, thresh = 0.05) {
  X <- close_comp(as.matrix(X)); n <- nrow(X); D <- ncol(X)
  b <- max(10L, floor(sub_ratio * n))
  freq <- array(0, dim = c(length(grid), D, D))
  for (rep in seq_len(B)) {
    Xs <- close_comp(X[sample.int(n, b), , drop = FALSE])
    for (g in seq_along(grid)) freq[g, , ] <- freq[g, , ] + estimator(Xs, grid[g])
  }
  freq <- freq / B
  inst <- apply(freq, 1, function(P) { xi <- 2 * P * (1 - P); mean(xi[upper.tri(xi)]) })
  inst_mono <- cummax(inst)
  ok <- which(inst_mono <= thresh)
  gstar <- if (length(ok)) max(ok) else 1L     # densest stable graph (smallest penalty under thresh)
  list(adj = estimator(X, grid[gstar]), penalty = grid[gstar], instability = inst, grid = grid)
}

#' Graphical-lasso edge estimator on a clr covariance (closure-friendly wrapper).
#'
#' Thin wrapper around \code{glasso_read_clr()} returning only the 0/1
#' adjacency, for use as an estimator closure in \code{stars_generic()}.
#'
#' @param S_clr D x D clr covariance (may be singular / rank D-1).
#' @param rho graphical-lasso penalty on the correlation scale.
#' @param ridge_frac ridge fraction for the singular-matrix convergence
#'   fallback (default 1e-3).
#' @return D x D 0/1 adjacency (zero diagonal).
#' @noRd
glasso_adj <- function(S_clr, rho, ridge_frac = 1e-3) {
  glasso_read_clr(S_clr, rho = rho, ridge_frac = ridge_frac)$adj
}

#' StARS penalty selection for the cellNetCoDa back ends.
#'
#' Selects the sparsity penalty by subsampling stability (StARS): for each
#' penalty the graph is refitted on \code{B} subsamples, the per-edge selection
#' frequency gives the edge instability \eqn{2\theta(1-\theta)} averaged over
#' edges, and the densest graph whose monotonised instability stays at or below
#' \code{thresh} is returned, refitted on the full data at the selected penalty.
#'
#' @param X compositional data (closed internally).
#' @param method back end, \code{"mb"} or \code{"glasso"}.
#' @param rho_grid glasso penalty grid (decreasing); required for \code{"glasso"}.
#' @param lambda_grid neighbourhood-selection penalty grid (decreasing);
#'   required for \code{"mb"}.
#' @param B number of subsamples (default 20).
#' @param sub_ratio subsample fraction (default 0.7).
#' @param thresh instability threshold (default 0.05).
#' @param k components for the cellwise detector (passed to \code{detect_cells}).
#' @param rule neighbourhood-combination rule for \code{"mb"}.
#' @param \dots forwarded to \code{detect_cells()}.
#' @return list with \code{adj}, \code{penalty}, \code{instability},
#'   \code{instability_mono}, \code{grid} and \code{method}.
#' @noRd
stars_select <- function(X, method = c("mb", "glasso"),
                         rho_grid = NULL, lambda_grid = NULL,
                         B = 20, sub_ratio = 0.7, thresh = 0.05, k = NULL, rule = "and", ...) {
  method <- match.arg(method)
  X <- close_comp(as.matrix(X)); n <- nrow(X); D <- ncol(X)
  grid <- if (method == "glasso") rho_grid else lambda_grid
  stopifnot(!is.null(grid))
  b <- max(10L, floor(sub_ratio * n))
  # accumulate edge-selection frequency per grid value across B subsamples
  freq <- array(0, dim = c(length(grid), D, D))
  for (rep in seq_len(B)) {
    idx <- sample.int(n, b)
    Xs  <- close_comp(X[idx, , drop = FALSE])
    det <- detect_cells(Xs, k = k, ...)
    rc  <- robust_ilr_cov(Xs, det)
    for (g in seq_along(grid)) {
      adj <- if (method == "glasso")
        precision_glasso(rc, rho = grid[g])$adj
      else
        precision_mb(Xs, det, rc, lambda = grid[g], rule = rule)$adj
      freq[g, , ] <- freq[g, , ] + adj
    }
  }
  freq <- freq / B
  # StARS edge instability: 2 * theta * (1 - theta), averaged over edges
  inst <- apply(freq, 1, function(P) {
    xi <- 2 * P * (1 - P); mean(xi[upper.tri(xi)])
  })
  # monotonise (StARS uses the running max as penalty decreases / density increases)
  inst_mono <- cummax(inst)
  ok <- which(inst_mono <= thresh)
  gstar <- if (length(ok)) max(ok) else 1L     # densest stable graph (standard StARS; smallest penalty under the instability threshold, grid ordered decreasing-penalty)
  # final fit on full data at the selected penalty
  detF <- detect_cells(X, k = k, ...); rcF <- robust_ilr_cov(X, detF)
  adjF <- if (method == "glasso") precision_glasso(rcF, rho = grid[gstar])$adj
          else precision_mb(X, detF, rcF, lambda = grid[gstar], rule = rule)$adj
  list(adj = adjF, penalty = grid[gstar], instability = inst, instability_mono = inst_mono,
       grid = grid, method = method)
}

## ---- detect_cells.R ----

#' Masking-resistant cellwise detection -> per-part soft weights (n x D).
#'
#' Flags cellwise outliers by a conditional-pivot test and grades them into
#' soft weights. For each part \eqn{j} the first ilr pivot coordinate isolating
#' that part,
#' \deqn{z^{(j)}_i = \sqrt{(D-1)/D}\,\bigl(\log x_{ij} - \tfrac{1}{D-1}\sum_{l\neq j}\log x_{il}\bigr),}
#' is robustly regressed on the (D-2) pivot coordinates of the remaining parts
#' (an MM regression via \code{robustbase::lmrob}, falling back to OLS), and a
#' cell is flagged when the robustly standardised residual exceeds \code{k_thr}.
#' The scale is re-estimated on the cells that survive the cutoff and iterated
#' (\code{iter} passes) to remove self-masking. Conditioning each part on the
#' rest explains away the network-driven covariance, so a cellwise perturbation
#' stands out against a smaller residual scale than a marginal clr detector
#' would see.
#'
#' The flagging cutoff \code{k_thr = 2.5} is the standard cellwise-detection
#' cutoff of cellMCD (Raymaekers and Rousseeuw, 2024). Detection (which cells)
#' is decoupled from grading (how hard to downweight them): clean cells keep
#' weight 1; flagged cells are graded by \code{"marginal"} (Tukey weight on the
#' standardised clr residual), \code{"projection"} (a correlation-adjusted
#' propagation-direction score mirroring \code{cellPcaCoDa}), or
#' \code{"isotropic"} (every flagged cell gets \code{w0}).
#'
#' The expensive part (the per-part \code{lmrob} regressions and the
#' \code{cellPcaCoDa} fit) depends only on the data and is factored into
#' \code{detect_cells_core()}; \code{grade_weights()} is the cheap
#' grading-to-weights step. Callers needing several gradings on the same data
#' should call \code{detect_cells_core()} once and \code{grade_weights()} per
#' grading.
#'
#' @param X     n x D composition matrix (positive entries; closed internally).
#' @param k     Number of PCA components passed to cellPcaCoDa (NULL = auto).
#' @param w0    Floor weight for flagged cells (default 0.05).
#' @param k_thr Conditional-residual flagging cutoff (default 2.5; cellMCD
#'   convention).
#' @param c_tukey  Tukey bisquare cutoff used only by the weight grading
#'   (default 3).
#' @param iter  Number of de-masking passes for the conditional scale (default
#'   3; >= 2 enables self-masking removal).
#' @param grading Weight-magnitude strategy for flagged cells: "marginal"
#'   (default), "projection" (correlation-adjusted propagation-direction score),
#'   or "isotropic" (all flagged cells get w0).
#' @param ...   Additional arguments forwarded to cellPcaCoDa.
#' @return list with \code{weights} (n x D soft weights), \code{flags} and the
#'   intermediate detection statistics produced by \code{detect_cells_core}.
#' @noRd
detect_cells <- function(X, k = NULL, w0 = 0.05, k_thr = 2.5, c_tukey = 3, iter = 3L,
                         grading = c("marginal", "projection", "isotropic"), ...) {
  grading <- match.arg(grading)
  core <- detect_cells_core(X, k = k, k_thr = k_thr, iter = iter, ...)
  W <- grade_weights(core, grading = grading, w0 = w0, c_tukey = c_tukey)
  list(weights = W, flags = core$flags, mu_clr = core$mu_clr,
       flags_cond = core$flags_cond, flags_pca = core$flags_pca,
       cond_abs = core$cond_abs, fit = core$fit)
}

#' Grading-independent core of \code{detect_cells} (the expensive part).
#'
#' Runs the masking-resistant conditional-pivot detector (per-part \code{lmrob}
#' regressions with iterative de-masking) and the \code{cellPcaCoDa} fit, takes
#' the union of their flags, and precomputes every residual statistic the three
#' gradings need: the marginal robustly standardised clr residual \code{Rstd}
#' and the correlation-adjusted propagation-direction score \code{proj_std}.
#' Depends only on \code{X} / \code{k} / \code{k_thr} / \code{iter}, not on
#' \code{grading}, \code{w0} or \code{c_tukey}, so its result can be cached once
#' per subsample.
#'
#' @param X     n x D composition matrix (positive entries; closed internally).
#' @param k     Number of PCA components passed to cellPcaCoDa (NULL = auto).
#' @param k_thr Conditional-residual flagging cutoff (default 2.5).
#' @param iter  Number of de-masking passes for the conditional scale (default 3).
#' @param ...   Additional arguments forwarded to cellPcaCoDa.
#' @return list with components used by \code{grade_weights}: \code{flags}
#'   (union, n x D logical), \code{mu_clr}, \code{flags_cond}, \code{flags_pca},
#'   \code{cond_abs}, \code{fit}, \code{Rstd} (n x D marginal standardised clr
#'   residual) and \code{proj_std} (n x D correlation-adjusted projection score).
#' @noRd
detect_cells_core <- function(X, k = NULL, k_thr = 2.5, iter = 3L, ...) {
  X <- close_comp(as.matrix(X)); n <- nrow(X); D <- ncol(X)
  if (D < 3) stop("detect_cells requires D >= 3 parts")
  has_rb <- requireNamespace("robustbase", quietly = TRUE)

  ## --- cellPcaCoDa: robust clr centre + its own (projection) flags (one detection source) ---
  fit <- cellPcaCoDa(as.data.frame(X), k = k, ...)
  flags_pca <- fit$cellflags                                # n x D logical
  mu_clr    <- fit$center_clr                               # length D robust clr centre

  ## --- Conditional-pivot detector (the masking-resistant source) ---
  L  <- log(X)                                              # n x D log-abundances
  Vr <- get_V(D - 1L)                                       # (D-1) x (D-2) basis for the "rest"
  cond_abs <- matrix(NA_real_, n, D)                        # |robust conditional std residual|
  for (j in seq_len(D)) {
    rest <- seq_len(D)[-j]
    # first pivot coordinate isolating part j vs the geometric mean of the rest
    z1 <- sqrt((D - 1) / D) * (L[, j] - rowMeans(L[, rest, drop = FALSE]))
    # (D-2) pivot coordinates among the rest (do not involve part j directly)
    Lr <- L[, rest, drop = FALSE]
    Zr <- (Lr - rowMeans(Lr)) %*% Vr                        # n x (D-2)
    fitj <- if (has_rb)
      tryCatch(suppressWarnings(robustbase::lmrob(z1 ~ Zr, setting = "KS2014")),
               error = function(e) NULL) else NULL
    r <- if (is.null(fitj)) stats::residuals(stats::lm(z1 ~ Zr)) else stats::residuals(fitj)
    s <- stats::mad(r, constant = 1.4826); if (s < 1e-8) s <- 1e-8
    rs <- (r - stats::median(r)) / s
    # iterative de-masking: re-estimate centre/scale on the cells that survive k_thr
    if (iter > 1L) for (it in seq_len(iter - 1L)) {
      keep <- abs(rs) <= k_thr
      if (sum(keep) > 10L) {
        s2 <- stats::mad(r[keep], constant = 1.4826); if (s2 < 1e-8) s2 <- 1e-8
        rs <- (r - stats::median(r[keep])) / s2
      }
    }
    cond_abs[, j] <- abs(rs)
  }
  flags_cond <- cond_abs > k_thr

  ## Union of detection sources (matches the original union design; cond is the workhorse).
  flags_union <- flags_cond | flags_pca

  ## --- clr residuals + robust per-column scale, for the WEIGHT GRADING only ---
  ## Both grading statistics (marginal Rstd, projection proj_std) are grading-independent and
  ## cheap, so we precompute BOTH here; grade_weights() just selects + applies the floor.
  R <- sweep(clr_mat(X), 2, mu_clr, "-")
  s_clr <- apply(R, 2, function(z) stats::mad(z, constant = 1.4826)); s_clr[s_clr < 1e-8] <- 1e-8
  Rstd  <- sweep(R, 2, s_clr, "/")                          # n x D (marginal grading input)

  ## projection grading input: correlation-adjusted propagation-direction score
  ## (mirrors cellPcaCoDa's detector; here applied to the robust-centred clr residuals).
  V  <- get_V(D)                                            # D x (D-1)
  Zr2 <- R %*% V                                            # n x (D-1) ilr residuals
  sz  <- apply(Zr2, 2, function(z) stats::mad(z, constant = 1.4826)); sz[sz < 1e-8] <- 1e-8
  Zs  <- sweep(Zr2, 2, sz, "/")
  proj_dirs <- matrix(0, nrow = D - 1L, ncol = D)
  for (j in seq_len(D)) {
    d_j <- as.numeric(crossprod(V, diag(D)[, j] - 1 / D))
    nrm <- sqrt(sum(d_j^2)); if (nrm < 1e-14) nrm <- 1
    proj_dirs[, j] <- d_j / nrm
  }
  proj_scores <- Zs %*% proj_dirs                           # n x D
  Rcor <- stats::cor(Zs); Rcor[!is.finite(Rcor)] <- 0; diag(Rcor) <- 1
  sigma_p <- sqrt(pmax(colSums((Rcor %*% proj_dirs) * proj_dirs), .Machine$double.eps))
  proj_std <- sweep(proj_scores, 2, sigma_p, "/")           # n x D (projection grading input)

  list(flags = flags_union, mu_clr = mu_clr,
       flags_cond = flags_cond, flags_pca = flags_pca, cond_abs = cond_abs, fit = fit,
       Rstd = Rstd, proj_std = proj_std)
}

#' Cheap grading-to-weights step: turn a \code{detect_cells_core} object into
#' n x D soft weights.
#'
#' Clean (unflagged) cells keep weight 1; flagged cells are downweighted by the
#' chosen grading and floored at \code{w0}. O(n*D); no \code{lmrob} /
#' \code{cellPcaCoDa}.
#'
#' @param core    Output of \code{detect_cells_core()}.
#' @param grading "marginal", "projection", or "isotropic".
#' @param w0      Floor weight for flagged cells (default 0.05).
#' @param c_tukey Tukey bisquare cutoff for the graded weight (default 3).
#' @return n x D weight matrix in [0, 1].
#' @noRd
grade_weights <- function(core, grading = c("marginal", "projection", "isotropic"),
                          w0 = 0.05, c_tukey = 3) {
  grading <- match.arg(grading)
  fl <- core$flags
  W <- matrix(1, nrow(fl), ncol(fl))
  if (grading == "marginal") {
    W[fl] <- tukey_w(core$Rstd[fl], c = c_tukey)
    W[fl] <- pmax(W[fl], w0)
  } else if (grading == "isotropic") {
    W[fl] <- w0
  } else {
    ## grading == "projection"
    W[fl] <- tukey_w(core$proj_std[fl], c = c_tukey)
    W[fl] <- pmax(W[fl], w0)
  }
  W
}
