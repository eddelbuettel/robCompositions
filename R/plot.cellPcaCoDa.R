#' Print method for cellwise-robust compositional PCA
#'
#' Prints a concise summary of a \code{cellPcaCoDa} object including the
#' number of retained components, convergence status, the number and
#' percentage of flagged cells, and the cumulative explained variability.
#'
#' @param x object of class \code{"cellPcaCoDa"} as returned by
#'   \code{\link{cellPcaCoDa}}
#' @param \dots additional arguments (currently ignored)
#'
#' @return Invisibly returns \code{x}.
#' @author Matthias Templ
#' @seealso \code{\link{cellPcaCoDa}}, \code{\link{summary.cellPcaCoDa}},
#'   \code{\link{plot.cellPcaCoDa}}
#' @keywords print
#' @export
#' @method print cellPcaCoDa
#' @examples
#' data(arcticLake)
#' res <- cellPcaCoDa(arcticLake, k = 2)
#' res
print.cellPcaCoDa <- function(x, ...) {
  cat("\nCellwise-robust PCA for compositional data\n")
  cat("-------------------------------------------\n")
  eV    <- x$eigenvalues / sum(x$eigenvalues_all)
  eVcum <- cumsum(eV)
  cat(sprintf("Components: %d\n", x$k))
  cat(sprintf("Converged:  %s (iterations: %d)\n", x$converged, x$iterations))
  cat(sprintf("Flagged cells: %d / %d (%.1f%%)\n",
              sum(x$cellflags), length(x$cellflags),
              100 * mean(x$cellflags)))
  cat("\nCumulative explained variability (retained components):\n")
  names(eVcum) <- paste0("Comp.", seq_along(eVcum))
  print(round(eVcum, 4))
  cat("\n")
  invisible(x)
}

#' Summary method for cellwise-robust compositional PCA
#'
#' Provides a more detailed summary than \code{\link{print.cellPcaCoDa}},
#' including individual and cumulative proportions of variance for all
#' eigenvalues as well as the contamination rate per variable.
#'
#' @param object object of class \code{"cellPcaCoDa"} as returned by
#'   \code{\link{cellPcaCoDa}}
#' @param \dots additional arguments (currently ignored)
#'
#' @return Invisibly returns a list with components \code{importance}
#'   (a matrix of variance proportions) and \code{flagrate} (per-variable
#'   contamination rate).
#' @author Matthias Templ
#' @seealso \code{\link{cellPcaCoDa}}, \code{\link{print.cellPcaCoDa}},
#'   \code{\link{plot.cellPcaCoDa}}
#' @keywords print
#' @export
#' @method summary cellPcaCoDa
#' @examples
#' data(arcticLake)
#' res <- cellPcaCoDa(arcticLake, k = 2)
#' summary(res)
summary.cellPcaCoDa <- function(object, ...) {
  stopifnot(inherits(object, "cellPcaCoDa"))

  ev_all <- object$eigenvalues_all
  prop   <- ev_all / sum(ev_all)
  cumprop <- cumsum(prop)
  importance <- rbind("Eigenvalue"         = ev_all,
                      "Proportion"         = prop,
                      "Cumulative Proportion" = cumprop)
  colnames(importance) <- paste0("Comp.", seq_along(ev_all))

  cat("\nCellwise-robust PCA for compositional data\n")
  cat("-------------------------------------------\n")
  cat(sprintf("Components retained: %d of %d\n", object$k, length(ev_all)))
  cat(sprintf("Converged: %s (%d iterations)\n", object$converged,
              object$iterations))
  cat(sprintf("Loss function: %s\n\n", object$rho))

  cat("Importance of components:\n")
  print(round(importance, 4))

  ## Per-variable flagging rate
  flagrate <- colMeans(object$cellflags)
  cat("\nPer-variable contamination rate:\n")
  print(round(flagrate, 4))
  cat("\n")

  invisible(list(importance = importance, flagrate = flagrate))
}

#' Plot method for cellwise-robust compositional PCA
#'
#' Provides diagnostic plots for a \code{cellPcaCoDa} object.
#' \code{which = 1} (default) draws a contamination heatmap of the raw-part
#' flags.
#' \code{which = 2} draws a compositional biplot in clr space.
#'
#' @param x object of class \code{"cellPcaCoDa"} as returned by
#'   \code{\link{cellPcaCoDa}}
#' @param y not used
#' @param \dots further arguments passed to \code{\link[graphics]{image}}
#'   (for \code{which = 1}) or \code{\link[stats]{biplot}} (for
#'   \code{which = 2})
#' @param which integer; \code{1} for contamination heatmap (default),
#'   \code{2} for compositional biplot
#'
#' @return Called for its side effect (a plot).  Returns \code{NULL}
#'   invisibly.
#' @author Matthias Templ
#' @seealso \code{\link{cellPcaCoDa}}, \code{\link{print.cellPcaCoDa}},
#'   \code{\link{summary.cellPcaCoDa}}
#' @keywords hplot
#' @importFrom graphics image axis box mtext
#' @importFrom stats biplot
#' @export
#' @method plot cellPcaCoDa
#' @examples
#' data(arcticLake)
#' set.seed(123)
#' x <- arcticLake
#' x[1, 2] <- x[1, 2] * 10
#' x <- constSum(x)
#' res <- cellPcaCoDa(x, k = 2)
#'
#' ## contamination heatmap
#' plot(res)
#'
#' ## biplot
#' plot(res, which = 2)
plot.cellPcaCoDa <- function(x, y, ..., which = 1) {
  if (which == 1) {
    ## --- Cellwise contamination heatmap ---
    flags <- x$cellflags * 1L
    nr <- nrow(flags)
    nc <- ncol(flags)

    ## image() expects a matrix with rows = x-axis and columns = y-axis,
    ## so we transpose and reverse the row order to get observations on
    ## the y-axis (top to bottom) and parts on the x-axis.
    img <- t(flags[nr:1, , drop = FALSE])

    ## axis tick positions (0--1 normalised)
    xat <- if (nc > 1) seq(0, 1, length.out = nc) else 0.5
    yat <- if (nr > 1) seq(0, 1, length.out = nr) else 0.5

    ## x-axis labels: part names or indices
    xlabs <- if (!is.null(colnames(x$cellflags))) {
      colnames(x$cellflags)
    } else {
      seq_len(nc)
    }
    ## y-axis labels: observation indices (top to bottom)
    ylabs_idx <- seq_len(nr)
    ## only show a subset if many observations (avoid crowding)
    show_y <- if (nr <= 40) seq_len(nr) else
      unique(c(1, round(seq(1, nr, length.out = 20)), nr))
    yat_sub <- yat[nr - show_y + 1]
    ylabs_sub <- show_y

    image(x = xat, y = yat, z = img,
          col = c("white", "red"),
          xlab = "", ylab = "",
          xaxt = "n", yaxt = "n", ...)
    axis(1, at = xat, labels = xlabs, las = 2, tick = FALSE)
    axis(2, at = yat_sub, labels = ylabs_sub, las = 1, tick = FALSE)
    mtext("Part", side = 1, line = 3)
    mtext("Observation", side = 2, line = 3)
    box()
    title(main = "Cellwise contamination flags")

  } else if (which == 2) {
    ## --- Biplot in clr space ---
    if (x$k < 2) {
      stop("biplot requires at least 2 components (k >= 2)")
    }
    biplot(x$scores[, 1:2], x$loadings[, 1:2],
           xlab = "Comp.1 (clr)", ylab = "Comp.2 (clr)",
           main = "Cellwise-robust compositional biplot", ...)
  } else {
    stop("'which' must be 1 or 2")
  }
  invisible(NULL)
}
