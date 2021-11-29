getCriterion <- function(mysample, dpnt, weight, phi) {
  D2dpnt <- sqrt((mysample$PC1 - dpnt$x1)^2 + (mysample$PC2 - dpnt$x2)^2)
  D <- as.matrix(dist(mysample[, c("x1", "x2")]))
  if (is.null(phi)) {
    diag(D) <- NA
    logdmin <- apply(D, MARGIN = 1, FUN = min, na.rm = TRUE) %>% log
    criterion_cur <- mean(-logdmin) + mean(D2dpnt) * weight
  } else {
    C <- variogramLine(vgm(model = "Exp", psill = 1, range = phi), dist_vector = D, covariance = TRUE)
    criterion_cur <- mean(C) + mean(D2dpnt) * weight
  }
}


permute <- function(mysample, candidates)  {
  unit.rand <- sample(nrow(mysample), size = 1)
  #remove selected unit from candidates
  candidates <- candidates[candidates$unit != mysample$unit[unit.rand], ]
  #replace selected sampling point by another randomly selected point from the same group
  units.dpnt <- which(candidates$dpnt == mysample$dpnt[unit.rand])
  candidates.dpnt <- candidates[units.dpnt, ]
  unit.candidate <- sample(nrow(candidates.dpnt), size = 1)
  mysample[unit.rand, ] <- candidates.dpnt[unit.candidate, ]
  mysample
}

# @ param mysample. Data frame of units of initial sample, with unit ID, 
#                   two spatial coordinates, two principal component scores, 
#                   and associated design-point 
# @ param candidates. Data frame with two or more candidates per design point.
#                     Data frame has same variables as mysample
# @ param dpnt. Data frame with coordinates of points of a rotatable central composite
#               design for two factors 
# @ param weight. Weight for the average distance of the sampling points to the
#                 associated design points (optional)  
# @ param phi. Distance parameter of exponential semivariogram (optional); 
#              if NULL (default value) the geometric criterion is used, 
#              else the model-based criterion is used
# @ param T_ini. Initial temperature of annealing schedule
# @ param coolingRate. Cooling rate of annealing schedule
# @ param maxPermuted. Multiplier for maximum number of permutations per chain
# @ param MaxNoChange. Stopping criterion. Maximum number of successive chains
#                      without change of the criterion
# @ return list with data frame of optimised sample, numeric with mimimised 
#          criterion, and numeric with trace of criterion 

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

  # set initial temperature
  T <- T_ini

  # compute the minimisation criterion MSSD for initial sample
  D2dpnt <- sqrt((mysample0$PC1 - dpnt$x1)^2 + (mysample0$PC2 - dpnt$x2)^2)
  D <- as.matrix(dist(mysample0[, c("x1", "x2")]))
  if (!is.null(phi)) {
    C <- variogramLine(vgm(model = "Exp", psill = 1, range = phi), dist_vector = D, covariance = TRUE)
    criterion_cur <- mean(C) + mean(D2dpnt) * weight
  } else {
    diag(D) <- NA
    logdmin <- apply(D, MARGIN = 1, FUN = min, na.rm = TRUE) %>% log(.)
    criterion_cur <- mean(-logdmin) + mean(D2dpnt) * weight
  }

  # Define structure for storing trace of criterion
  trace <- NULL

  # initialize number of zero changes of objective function
  nNoChange <- 0

  mysample_cur <- mysample0

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
