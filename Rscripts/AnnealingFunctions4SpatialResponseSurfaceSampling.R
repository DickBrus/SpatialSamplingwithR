#' Validate Column Names
#'
#' This function checks whether all columns of data frame with sampling points
#' and of data frame with candidate sampling points are available
#'
#' @param x data.frame to validate
#'
#' @return validated data frame
#'
#' @importFrom dplyr "%>%"
#' @importFrom stringr str_glue
validate_colnames <- function(x) {

  required_column_names <- c("point_id", "x", "y", "PC1", "PC2", "dpnt")
  missing_column_names <- required_column_names %>%
    setdiff(names(x))
  n_missing_column_names <- length(missing_column_names)
  if (n_missing_column_names == 1L) {
    stop(str_glue("Column {sQuote(missing_column_names)} is missing"), call. = FALSE)
  }
  if (n_missing_column_names  > 1L) {
    stop(paste0("Columns ", toString(sQuote(missing_column_names)), " are missing"), call. = FALSE)
  }

  return(x)
}

#' Minimisation Criterion
#'
#' This function computes the minimisation criterion
#' (value of objective function).
#'
#' @param mysample \code{\link{data.frame}} of units of proposed sample.
#' @param dpnt \code{\link{data.frame}} with coordinates of design points of a
#'   rotatable central composite design for two factors.
#' @param weight weight for the average distance of the sampling points to the
#'   associated design points (optional).
#' @param phi distance parameter of exponential semivariogram (optional);
#'   if \code{NULL} (default value) the geometric criterion is used,
#'   else the model-based criterion is used.
#' @return minimisation criterion value of proposed sample.
#'
#' @importFrom stats dist
#' @importFrom gstat variogramLine vgm
getCriterion <- function(mysample, dpnt, weight, phi = NULL) {
  D2dpnt <- sqrt((mysample$PC1 - dpnt$x1)^2 + (mysample$PC2 - dpnt$x2)^2)
  D <- as.matrix(dist(mysample[, c("x", "y")]))
  if (is.null(phi)) {
    diag(D) <- NA
    dmin <- apply(D, MARGIN = 1, FUN = min, na.rm = TRUE)
    logdmin <- log(dmin)
    criterion_cur <- mean(-logdmin) + mean(D2dpnt) * weight
  } else {
    C <- variogramLine(vgm(model = "Exp", psill = 1, range = phi),
                       dist_vector = D, covariance = TRUE)
    criterion_cur <- mean(C) + mean(D2dpnt) * weight
  }
}



#' Permutation
#'
#' This function permutes a current sample by replacing one randomly selected
#' unit by a unit randomly selected from the remaining candidate units of the
#' same design point.
#'
#' @param mysample \code{\link{data.frame}} of units of current sample.
#' @param candidates \code{\link{data.frame}} with candidate units per design.
#'   point. This \code{\link{data.frame}} has same variables as \code{mysample}.
#'
#' @return \code{\link{data.frame}} with units of permuted sample.
permute <- function(mysample, candidates)  {
  unit_rand <- sample(nrow(mysample), size = 1)
  #remove selected unit from candidates
  candidates <- candidates[candidates$point_id != mysample$point_id[unit_rand], ]
  #replace selected sampling point by another randomly selected point from the same group
  units_dpnt <- which(candidates$dpnt == mysample$dpnt[unit_rand])
  cnd_dpnt <- candidates[units_dpnt, ]
  unit_cnd <- sample(nrow(cnd_dpnt), size = 1)
  mysample[unit_rand, ] <- cnd_dpnt[unit_cnd, ]
  mysample
}



#' Simulated Annealing
#'
#' This function optimises the sampling pattern of a
#' central composite response surface design sample.
#'
#' @param mysample \code{\link{data.frame}} of units of initial sample,
#'   with point_id, two spatial coordinates, two principal component scores,
#'   and associated design-point.
#' @param candidates \code{\link{data.frame}} with candidate units per design
#'   point. \code{\link{data.frame}} has same variables as \code{mysample}.
#' @param dpnt \code{\link{data.frame}} with coordinates of design points of a
#'   rotatable central composite design for two factors.
#' @param weight weight for the average distance of the sampling points to the
#    associated design points (optional).
#' @param phi distance parameter of exponential semivariogram (optional);
#'   if \code{NULL} (default value) the geometric criterion is used,
#'   else the model-based criterion is used.
#' @param T_ini initial temperature of annealing schedule.
#' @param coolingRate cooling rate of annealing schedule.
#' @param maxPermuted maximum number of permutations per chain.
#' @param maxNoChange stopping criterion. Maximum number of successive chains
#'   without change of the criterion.
#' @param verbose \code{\link{logical}} for printing the temperature, the
#'   criterion, the number of accepted proposals per chain,
#'   and the number of improvements per chain during the optimisation.
#'
#' @return \code{\link{list}} with \code{\link{data.frame}} of optimised sample,
#'    \code{\link{numeric}} with minimised criterion, and \code{\link{numeric}}
#'     with trace of criterion.
#'
#' @importFrom stats dist runif
#'
#' @export
anneal <- function(mysample,
                 candidates,
                 dpnt,
                 weight = 0,
                 phi = NULL,
                 T_ini  =  1,
                 coolingRate  =  0.95,
                 maxPermuted = 20 * nrow(mysample),
                 maxNoChange = 20,
                 verbose = getOption("verbose")) {
  #check columns of mysample
  mysample <- validate_colnames(mysample)


  #check columns of candidates
  candidates <- validate_colnames(candidates)

  #check columns of dpnt
  if (!("x" %in% names(dpnt))) {
    stop("x must be in data frame dpnt", call. = FALSE)
  }
  if (!("y" %in% names(dpnt))) {
    stop("y must be in data frame dpnt", call. = FALSE)
  }

  # set initial temperature
  T <- T_ini

  # compute the minimisation criterion MSSD for initial sample
  D2dpnt <- sqrt((mysample$PC1 - dpnt$x1)^2 + (mysample$PC2 - dpnt$x2)^2)
  D <- as.matrix(dist(mysample[, c("x", "y")]))
  if (is.null(phi)) {
    diag(D) <- NA
    dmin <- apply(D, MARGIN = 1, FUN = min, na.rm = TRUE)
    logdmin <- log(dmin)
    criterion_cur <- mean(-logdmin) + mean(D2dpnt) * weight
  } else {
    C <- variogramLine(vgm(model = "Exp", psill = 1, range = phi), dist_vector = D, covariance = TRUE)
    criterion_cur <- mean(C) + mean(D2dpnt) * weight
  }

  # Define structure for storing trace of criterion
  trace <- NULL

  # initialize number of zero changes of objective function
  nNoChange <- 0

  mysample_cur <- mysample

  # start cooling loop
  repeat{

    # initialize number of accepted configurations
    nAccepted <- 0

    # initialize number of permuted configurations
    nPermuted <- 0

    # initialize number of improved configurations
    nImproved <- 0

    # start permutation loop
    repeat {

      # increase the number of permutations
      nPermuted <- nPermuted + 1

      # propose new stratum boundaries
      mysample_p <- permute(mysample_cur, candidates)

      # compute the criterion of this new stratification
      criterion_p <- getCriterion(mysample_p, dpnt, weight, phi)

      # accept/reject proposal by means of Metropolis criterion
      delta <- criterion_p - criterion_cur
      if (delta < 0) {
        nImproved <- nImproved + 1
        p <- 1 # always accept improvements
      } else {
        p <- exp(-delta / T) # use Boltzmann to judge if deteriorations should be accepted
      }
      u <- runif(n = 1) # draw uniform deviate
      if (u < p) { # accept proposal
        nAccepted <- nAccepted + 1
        mysample_cur <- mysample_p
        criterion_cur <- criterion_p
      }
      # are conditions met to lower temperature?
      lowerTemperature <- (nPermuted == maxPermuted)
      if (lowerTemperature) {
        if (nImproved == 0) {
          nNoChange <- nNoChange + 1
          }
        else {
          nNoChange <- 0
          }
        trace <- rbind(trace, criterion_cur)
        break
      }
    }

    if (verbose) {
      cat(
        format(Sys.time()), "|",
        sprintf("T = %e  C = %e  permuted = %d  accepted = %d  improved = %d  acceptance rate = %f  \n",
                T, criterion_cur, nPermuted, nAccepted, nImproved, nAccepted / nPermuted)
      )
    }

    # check on convergence
    if (nNoChange == maxNoChange) {
      break
    }

    # lower temperature
    T <- coolingRate * T
  }

  # return result
  list(
    mysample = mysample_cur, criterion = criterion_cur, trace = trace)
}
