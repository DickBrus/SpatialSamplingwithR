#' Local Mean Variance Neighbors and Weights
#'
#' This function calculates the index values of neighboring points and
#' associated weights required by the local mean variance estimator.
#'
#' @param x \code{\link{vector}} of x-coordinates for location of the sample
#'   points.
#' @param y \code{\link{vector}} of y-coordinates for location of the sample
#'   points.
#' @param prb \code{\link{vector}} of inclusion probabilities for the sample
#'   points.
#' @param nbh Number of neighboring points to use in the calculations.
#'
#' @return If \code{\link{ginv}} fails to return valid output,
#'   a \code{\link{NULL}} object.  Otherwise, a \code{\link{list}}
#'   containing two elements: a \code{\link{matrix}} named \code{ij} composed of
#'   the index values of neighboring points and a vector named \code{gwt}
#'   composed of weights.
#'
#' @importFrom MASS ginv
#'
#' @author  Tom Kincaid \email{Kincaid.Tom@@epa.gov}
#'
#' @export
localmean_weight <- function(x, y, prb, nbh = 4) {

  # Assign tne number of points

  n <- length(x)

  # Calculate indices of nearest neighbors

  dst <- as.matrix(dist(cbind(x, y), diag = TRUE, upper = TRUE))
  idx <- apply(dst, 2, order)[1:nbh, ]

  # Make neighbors symmetric

  jdx <- rep(1:n, rep(nbh, n))
  kdx <- unique(c((jdx - 1) * n + idx, (idx - 1) * n + jdx)) - 1
  ij <- cbind((kdx) %/% n + 1, (kdx) %% n + 1)
  ij <- ij[order(ij[, 1], dst[ij]), ]

  # Apply linear taper to the  inverse probability weights

  gct <- tabulate(ij[, 1])
  gwt <- numeric(0)
  for (i in 1:n) {
    gwt <- c(gwt, 1 - (1:gct[i] - 1) / (gct[i]))
  }
  gwt <- gwt / prb[ij[, 2]]

  # Normalize to make true average

  smwt <- sapply(split(gwt, ij[, 1]), sum)
  gwt <- gwt / smwt[ij[, 1]]
  smwt <- sapply(split(gwt, ij[, 2]), sum)

  # Make weights doubly stochastic

  hij <- matrix(0, n, n)
  hij[ij] <- 0.5
  a22 <- try(ginv(diag(gct / 2) - hij %*% diag(2 / gct) %*% hij), TRUE)
  if ("try-error" %in% class(a22)) {
    return(NULL)
  }
  a21 <- -diag(2 / gct) %*% hij %*% a22
  lm <- a21 %*% (1 - smwt)
  gm <- a22 %*% (1 - smwt)
  gwt <- (lm[ij[, 1]] + gm[ij[, 2]]) / 2 + gwt

  # Return the results

  list(ij = ij, gwt = gwt)
}
