###############################################################################
# Function: localmean_var (not exported)
# Programmers: Don Stevens and Tom Kincaid
# Date: October 17, 2000
#
#' Internal Function: Local Mean Variance Estimator
#'
#' This function calculates the local mean variance estimator.
#'
#' @param z Vector of weighted response values or weighted residual values for
#'  the sample points.
#'
#' @param weight_1st List from the local mean weight function containing two
#'  elements: a matrix named \code{ij} composed of the index values of neighboring
#'  points and a vector named \code{gwt} composed of weights.
#'
#' @return The local mean estimator of the variance.
#'
#' @author Tom Kincaid \email{Kincaid.Tom@@epa.gov}
#'
#' @keywords survey
#'
#' @noRd
###############################################################################

localmean_var <- function(z, weight_1st) {

  # Calculate local means

  zb <- sapply(split(z[weight_1st$ij[, 2]] * weight_1st$gwt, weight_1st$ij[, 1]), sum)

  # Calculate the variance estimate

  lmvar <- sum(weight_1st$gwt * (z[weight_1st$ij[, 2]] - zb[weight_1st$ij[, 1]])^2)

  # Return the variance estimate

  lmvar
}
