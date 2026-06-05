#' Print method for cellwise-robust compositional networks
#'
#' Prints a concise summary of a \code{cellNetCoDa} object: the back end used,
#' the number of nodes (parts) and selected edges, the graph density, the
#' StARS-selected penalty, and the share of flagged cells.
#'
#' @param x object of class \code{"cellNetCoDa"} as returned by
#'   \code{\link{cellNetCoDa}}
#' @param \dots additional arguments (currently ignored)
#'
#' @return Invisibly returns \code{x}.
#' @author Matthias Templ
#' @seealso \code{\link{cellNetCoDa}}, \code{\link{summary.cellNetCoDa}},
#'   \code{\link{plot.cellNetCoDa}}
#' @keywords print
#' @export
#' @method print cellNetCoDa
#' @examples
#' \donttest{
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   data(expenditures)
#'   set.seed(1)
#'   net <- cellNetCoDa(expenditures, method = "mb",
#'                      lambda_grid = c(0.4, 0.2, 0.1, 0.05), B = 10, k = 2)
#'   net
#' }
#' }
print.cellNetCoDa <- function(x, ...) {
  D    <- nrow(x$adj)
  poss <- D * (D - 1) / 2
  back <- switch(x$method,
                 mb     = "Meinshausen-Buhlmann (neighbourhood selection)",
                 glasso = "graphical lasso (clr covariance)",
                 x$method)
  cat("\nCellwise-robust network for compositional data\n")
  cat("----------------------------------------------\n")
  cat(sprintf("Back end:         %s\n", back))
  cat(sprintf("Nodes (parts):    %d\n", D))
  cat(sprintf("Edges:            %d of %d (density %.1f%%)\n",
              x$n_edges, poss, 100 * x$n_edges / poss))
  cat(sprintf("Selected penalty: %.4g\n", x$penalty))
  if (!is.null(x$flags)) {
    cat(sprintf("Flagged cells:    %d / %d (%.1f%%)\n",
                sum(x$flags), length(x$flags), 100 * mean(x$flags)))
  }
  cat("\n")
  invisible(x)
}

#' Summary method for cellwise-robust compositional networks
#'
#' Provides a more detailed summary than \code{\link{print.cellNetCoDa}}: the
#' per-node degree of the estimated graph, the per-part cellwise contamination
#' rate, and the list of selected edges.
#'
#' @param object object of class \code{"cellNetCoDa"} as returned by
#'   \code{\link{cellNetCoDa}}
#' @param \dots additional arguments (currently ignored)
#'
#' @return Invisibly returns a list with components \code{degree} (incident
#'   edges per node), \code{density}, \code{flagrate} (per-part contamination
#'   rate), \code{edges} (a data.frame of selected edges), \code{penalty} and
#'   \code{method}.
#' @author Matthias Templ
#' @seealso \code{\link{cellNetCoDa}}, \code{\link{print.cellNetCoDa}},
#'   \code{\link{plot.cellNetCoDa}}
#' @keywords print
#' @export
#' @method summary cellNetCoDa
#' @examples
#' \donttest{
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   data(expenditures)
#'   set.seed(1)
#'   net <- cellNetCoDa(expenditures, method = "mb",
#'                      lambda_grid = c(0.4, 0.2, 0.1, 0.05), B = 10, k = 2)
#'   summary(net)
#' }
#' }
summary.cellNetCoDa <- function(object, ...) {
  stopifnot(inherits(object, "cellNetCoDa"))
  adj  <- object$adj
  D    <- nrow(adj)
  poss <- D * (D - 1) / 2
  labs <- if (!is.null(object$taxa)) object$taxa
          else if (!is.null(rownames(adj))) rownames(adj)
          else paste0("T", seq_len(D))

  deg  <- rowSums(adj); names(deg) <- labs
  dens <- object$n_edges / poss

  fl <- object$flags
  flagrate <- if (!is.null(fl)) {
    if (is.null(colnames(fl)) && ncol(fl) == D) colnames(fl) <- labs
    colMeans(fl)
  } else NULL

  ut <- which(upper.tri(adj) & adj != 0, arr.ind = TRUE)
  edges <- if (nrow(ut))
    data.frame(from = labs[ut[, 1]], to = labs[ut[, 2]],
               stringsAsFactors = FALSE)
  else data.frame(from = character(0), to = character(0))

  back <- switch(object$method,
                 mb     = "Meinshausen-Buhlmann (neighbourhood selection)",
                 glasso = "graphical lasso (clr covariance)",
                 object$method)

  cat("\nCellwise-robust network for compositional data\n")
  cat("----------------------------------------------\n")
  cat(sprintf("Back end:         %s\n", back))
  cat(sprintf("Nodes (parts):    %d\n", D))
  cat(sprintf("Edges:            %d of %d (density %.1f%%)\n",
              object$n_edges, poss, 100 * dens))
  cat(sprintf("Selected penalty: %.4g\n\n", object$penalty))

  cat("Node degree (incident edges), decreasing:\n")
  print(sort(deg, decreasing = TRUE))

  if (!is.null(flagrate)) {
    cat("\nPer-part contamination rate:\n")
    print(round(flagrate, 4))
  }

  if (nrow(edges)) {
    if (nrow(edges) <= 20) {
      cat("\nEdges:\n"); print(edges, row.names = FALSE)
    } else {
      cat(sprintf("\n%d edges (see the $edges element of the returned summary).\n",
                  nrow(edges)))
    }
  } else cat("\nNo edges selected.\n")
  cat("\n")

  invisible(list(degree = deg, density = dens, flagrate = flagrate,
                 edges = edges, penalty = object$penalty,
                 method = object$method))
}

#' Plot method for cellwise-robust compositional networks
#'
#' Diagnostic plots for a \code{cellNetCoDa} object.
#' \code{which = 1} (default) draws the estimated network with the parts as
#' nodes on a circular layout (node size proportional to degree).
#' \code{which = 2} draws the cellwise contamination heatmap of the detector
#' flags (parts on the x-axis, observations on the y-axis).
#' \code{which = 3} draws the StARS instability path used for penalty
#' selection (raw and monotonised instability over the penalty grid, with the
#' selection threshold and the chosen penalty marked).
#'
#' @param x object of class \code{"cellNetCoDa"} as returned by
#'   \code{\link{cellNetCoDa}}
#' @param y not used
#' @param \dots further arguments passed to \code{\link[graphics]{segments}}
#'   (for \code{which = 1}) or \code{\link[graphics]{image}} (for
#'   \code{which = 2}); ignored for \code{which = 3}
#' @param which integer; \code{1} for the network graph (default), \code{2}
#'   for the contamination heatmap, \code{3} for the StARS instability path
#'
#' @return Called for its side effect (a plot).  Returns \code{NULL}
#'   invisibly.
#' @author Matthias Templ
#' @seealso \code{\link{cellNetCoDa}}, \code{\link{print.cellNetCoDa}},
#'   \code{\link{summary.cellNetCoDa}}
#' @keywords hplot
#' @importFrom graphics plot.new plot.window segments points lines text title
#'   par image axis box mtext abline
#' @export
#' @method plot cellNetCoDa
#' @examples
#' \donttest{
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   data(expenditures)
#'   set.seed(1)
#'   net <- cellNetCoDa(expenditures, method = "mb",
#'                      lambda_grid = c(0.4, 0.2, 0.1, 0.05), B = 10, k = 2)
#'   plot(net)                # network graph
#'   plot(net, which = 2)     # contamination heatmap
#'   plot(net, which = 3)     # StARS instability path
#' }
#' }
plot.cellNetCoDa <- function(x, y, ..., which = 1) {
  if (which == 1) {
    ## --- Circular network plot of the estimated graph ---
    adj <- x$adj; D <- nrow(adj)
    labs <- if (!is.null(x$taxa)) x$taxa
            else if (!is.null(rownames(adj))) rownames(adj)
            else paste0("T", seq_len(D))
    ang <- pi / 2 - 2 * pi * (seq_len(D) - 1) / D   # node 1 at top, clockwise
    px  <- cos(ang); py <- sin(ang)

    op <- par(mar = c(1, 1, 3, 1)); on.exit(par(op))
    plot.new()
    plot.window(xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4), asp = 1)
    title(main = sprintf("cellNetCoDa graph (%s, %d edges)",
                         x$method, x$n_edges))

    ut <- which(upper.tri(adj) & adj != 0, arr.ind = TRUE)
    if (nrow(ut))
      segments(px[ut[, 1]], py[ut[, 1]], px[ut[, 2]], py[ut[, 2]],
               col = "grey40", lwd = 1.5, ...)

    deg  <- rowSums(adj)
    cexn <- 2 + 2.5 * deg / max(1, max(deg))         # node size ~ degree
    points(px, py, pch = 21, bg = "steelblue", col = "white", cex = cexn)
    text(1.2 * px, 1.2 * py, labels = labs, cex = 0.8, xpd = NA)

  } else if (which == 2) {
    ## --- Cellwise contamination heatmap (parts x observations) ---
    flags <- x$flags * 1L
    nr <- nrow(flags); nc <- ncol(flags)
    img <- t(flags[nr:1, , drop = FALSE])
    xat <- if (nc > 1) seq(0, 1, length.out = nc) else 0.5
    yat <- if (nr > 1) seq(0, 1, length.out = nr) else 0.5
    xlabs <- if (!is.null(colnames(flags))) colnames(flags)
             else if (!is.null(x$taxa)) x$taxa else seq_len(nc)
    show_y <- if (nr <= 40) seq_len(nr) else
              unique(c(1, round(seq(1, nr, length.out = 20)), nr))
    yat_sub <- yat[nr - show_y + 1]

    image(x = xat, y = yat, z = img, col = c("white", "red"),
          xlab = "", ylab = "", xaxt = "n", yaxt = "n", ...)
    axis(1, at = xat, labels = xlabs, las = 2, tick = FALSE)
    axis(2, at = yat_sub, labels = show_y, las = 1, tick = FALSE)
    mtext("Part", side = 1, line = 3)
    mtext("Observation", side = 2, line = 3)
    box()
    title(main = "Cellwise contamination flags")

  } else if (which == 3) {
    ## --- StARS instability path ---
    inst <- x$instability
    if (is.null(inst)) stop("no StARS instability stored in this object")
    has_grid  <- !is.null(x$grid) && length(x$grid) == length(inst)
    xv        <- if (has_grid) x$grid else seq_along(inst)
    xlb       <- if (has_grid) "penalty" else "grid index"
    inst_mono <- cummax(inst)

    op <- par(mar = c(4.5, 4.5, 3, 1)); on.exit(par(op))
    plot.new()
    ## penalties run high -> low (graph sparse -> dense): show that left -> right
    xr <- if (has_grid) rev(range(xv)) else range(xv)
    plot.window(xlim = xr,
                ylim = range(c(0, inst, inst_mono, x$thresh), na.rm = TRUE))
    segments(xv, 0, xv, inst, col = "grey80")
    points(xv, inst, pch = 16, col = "steelblue")
    lines(xv, inst_mono, col = "steelblue", lwd = 2)
    if (!is.null(x$thresh)) abline(h = x$thresh, lty = 2, col = "grey40")
    if (has_grid && !is.null(x$penalty)) abline(v = x$penalty, lty = 3, col = "red")
    axis(1); axis(2, las = 1); box()
    title(main = "StARS instability path", xlab = xlb,
          ylab = "edge instability")

  } else {
    stop("'which' must be 1 (network), 2 (contamination heatmap) or 3 (StARS instability path)")
  }
  invisible(NULL)
}
