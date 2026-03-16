#' Generate cellwise-contaminated compositional data
#'
#' Introduces cellwise contamination on the simplex by replacing individual
#' raw parts with contaminated values and re-closing the composition.
#'
#' For each cell independently, contamination occurs with probability
#' \code{epsilon}. Three contamination types are available:
#' \code{"multiplicative"} multiplies the cell by \code{delta},
#' \code{"replacement"} sets the cell to \code{delta}, and
#' \code{"additive"} adds \code{delta} times the row median.
#' After contamination, affected rows are re-closed to preserve the original
#' row sums.
#'
#' @param x clean compositional data (n x D matrix or data.frame with positive values)
#' @param epsilon per-cell contamination probability (default: 0.05)
#' @param delta contamination factor (default: 10)
#' @param type \code{"multiplicative"} (x_j * delta), \code{"replacement"}
#'   (x_j = delta), or \code{"additive"} (x_j + delta * median(row))
#' @param seed random seed for reproducibility
#'
#' @return A list with components:
#'   \item{x_contaminated}{the contaminated compositions (re-closed to same row sums)}
#'   \item{flags}{n x D logical matrix indicating which cells were contaminated}
#'   \item{epsilon}{contamination rate used}
#'
#' @author Matthias Templ
#' @keywords manip
#' @export
#' @examples
#' data(arcticLake)
#' cont <- contaminate_simplex(arcticLake, epsilon = 0.1, delta = 10, seed = 42)
#' sum(cont$flags)  # number of contaminated cells
#' rowSums(cont$x_contaminated)  # should match rowSums(arcticLake)
contaminate_simplex <- function(x, epsilon = 0.05, delta = 10,
                                type = "multiplicative", seed = NULL) {
  ## input validation
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data.frame")
  }
  x <- as.matrix(x)
  if (!is.numeric(x) || any(x <= 0, na.rm = TRUE)) {
    stop("'x' must contain only positive numeric values")
  }
  if (epsilon < 0 || epsilon > 1) {
    stop("'epsilon' must be between 0 and 1")
  }
  type <- match.arg(type, choices = c("multiplicative", "replacement", "additive"))

  n <- nrow(x)
  D <- ncol(x)

  ## set seed if provided
  if (!is.null(seed)) set.seed(seed)

  ## generate contamination flags
  flags <- matrix(runif(n * D) < epsilon, nrow = n, ncol = D)

  ## early return when no cells are flagged
  if (!any(flags)) {
    return(list(x_contaminated = x, flags = flags, epsilon = epsilon))
  }

  ## store original row sums for re-closure
  orig_rowsums <- rowSums(x)

  ## apply contamination
  x_cont <- x
  if (type == "multiplicative") {
    x_cont[flags] <- x[flags] * delta
  } else if (type == "replacement") {
    x_cont[flags] <- delta
  } else if (type == "additive") {
    row_medians <- apply(x, 1, median)
    med_mat <- matrix(row_medians, nrow = n, ncol = D)
    x_cont[flags] <- x[flags] + delta * med_mat[flags]
  }

  ## re-close only rows that have at least one contaminated cell
  affected <- which(rowSums(flags) > 0)
  if (length(affected) > 0) {
    x_cont[affected, ] <- constSum(x_cont[affected, , drop = FALSE],
                                   const = 1) *
      orig_rowsums[affected]
  }

  ## preserve column names
  colnames(x_cont) <- colnames(x)

  return(list(x_contaminated = x_cont, flags = flags, epsilon = epsilon))
}
