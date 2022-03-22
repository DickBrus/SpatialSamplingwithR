#' Validate input files
#'
#' This function checks various input files
#'
#' @param esample \code{\link[sp]{SpatialPoints}} used for prediction at the
#'   evaluation points or used for evaluation.
#' @param model \pkg{gstat} model type of a priori semivariogram model.
#' @param thetas parameters of semivariogram model.
#' @param perturbation proportion of a semivariogram parameter value added to
#'   that value.
validate <- function(esample, model, thetas, perturbation) {
  if (!(class(esample) %in% c("SpatialPoints"))) {
    stop("esample must be of class SpatialPoints", call. = FALSE)
  }
  if (!(model %in% c("Exp", "Sph"))) {
    stop("model must be a character: Exp or Sph", call. = FALSE)
  }
  if (length(thetas) > 2) {
    stop("numeric thetas must be of length 2", call. = FALSE)
  }
  if (thetas[1] > 1) {
    stop("first value of thetas (ratio of spatial dependence) must be < 1", call. = FALSE)
  }
  if (perturbation > 1) {
    stop("perturbation must be < 1, say < 0.05", call. = FALSE)
  }
}

#' Fisher Information Matrix
#'
#' This function computes the Fisher information matrix.
#'
#' @param A correlation matrix of sampling points.
#' @param D distance matrix of sampling points.
#' @param model \pkg{gstat} model type of a priori semivariogram model.
#' @param thetas parameters of semivariogram model.
#' @param perturbation proportion of a semivariogram parameter value added to
#'   that value.
#'
#' @return Fisher information matrix.
#'
#' @seealso \code{\link[gstat]{vgm}}
#'
#' @importFrom sp "coordinates<-" spDists
#' @importFrom gstat variogramLine vgm
FIM <- function(A, D, model, thetas, perturbation) {
  pA <- dA <- list()
  for (i in seq_len(length(thetas))) {
    thetas.pert <- thetas
    thetas.pert[i] <- (1 + perturbation) * thetas[i]
    vgmodel.pert <- vgm(model = model, psill = thetas.pert[1], range = thetas.pert[2], nugget = 1 - thetas.pert[1])
    pA[[i]] <- variogramLine(vgmodel.pert, dist_vector = D, covariance = TRUE)
    dA[[i]] <- (pA[[i]] - A) / (thetas[i] * perturbation)
  }
  cholA <- try(chol(A), silent = TRUE)
  if (inherits(cholA, "try-error")) {
    return(Inf)
  } else {
    I <- matrix(0, length(thetas), length(thetas))
    for (i in seq_len(length(thetas))) {
      m_i <- solve(A, dA[[i]])
      for (j in i:length(thetas)) {
        m_j <- solve(A, dA[[j]])
        I[i, j] <- I[j, i] <- 0.5 * sum(diag(m_i %*% m_j))
      }
    }
  }
  return(I)
}

#' Determinant Information Matrix
#'
#' This function computes the log of the reciprocal determinant of Fisher
#' information matrix.
#'
#' @param points \code{\link{data.frame}} or \code{SpatialPoints(DataFrame)}
#'   with coordinates of sampling points.
#' @param model \pkg{gstat} model type of a priori semivariogram model.
#' @param thetas parameters of semivariogram model.
#' @param perturbation proportion of a semivariogram parameter value added to
#'   that value.
#'
#' @return log of the reciprocal determinant of Fisher information matrix.
#'
#' @seealso \code{\link[gstat]{vgm}}
#'
#' @importFrom sp "coordinates<-" spDists
#' @importFrom gstat variogramLine vgm
#'
#' @export
logdet <- function(points, model, thetas, perturbation = 0.01)  {
  if (!(model %in% c("Exp", "Sph"))) {
    stop("model must be a character: Exp or Sph", call. = FALSE)
  }
  if (length(thetas) > 2) {
    stop("numeric thetas must be of length 2", call. = FALSE)
  }
  if (thetas[1] > 1) {
    stop("first value of thetas (ratio of spatial dependence) must be < 1", call. = FALSE)
  }
  if (perturbation > 1) {
    stop("perturbation must be < 1, say < 0.05", call. = FALSE)
  }

  D <- spDists(points)
  vgmodel <- vgm(model = model, psill = thetas[1], range = thetas[2], nugget = 1 - thetas[1])
  A <- variogramLine(vgmodel, dist_vector = D, covariance = TRUE)
  cholA <- try(chol(A), silent = TRUE)
  if (inherits(cholA, "try-error")) {
    return(Inf)
  } else {
    I <- FIM(A, D, model, thetas, perturbation)
    cholI <- try(chol(I), silent = TRUE)
    if (inherits(cholI, "try-error")) {
      return(Inf)
    } else {
      logdet <- -determinant(I, logarithm = TRUE)$modulus
      logdet
    }
  }
}


#' Mean Variance of Kriging Variance
#'
#' @param points \code{\link{data.frame}} or \code{SpatialPoints(DataFrame)}
#'   with coordinates of sampling points.
#' @param psample \code{\link[sp]{SpatialPoints}} used for prediction at the
#'   evaluation points.
#' @param esample \code{\link[sp]{SpatialPoints}} used for evaluation.
#' @param model \pkg{gstat} model type of a priori semivariogram model.
#' @param thetas parameters of semivariogram model.
#' @param perturbation proportion of a semivariogram parameter value added to
#'   that value.
#' @return mean variance of the kriging variance.
#'
#' @importFrom sp "coordinates<-" spDists
#' @importFrom gstat variogramLine vgm
#'
#' @export
MVKV <- function(points, psample, esample, model, thetas, perturbation = 0.01) {
  validate(esample = esample, model, thetas, perturbation)
  if (!(class(psample) %in% c("SpatialPoints"))) {
    stop("psample must be of class SpatialPoints", call. = FALSE)
  }
  points <- as.data.frame(points)
  nobs <- nrow(points)
  coordinates(points) <- ~x + y

  #compute distance matrix and correlation matrix of sampling points used for estimating the semivariogram
  D <- spDists(points)
  vgmodel <- vgm(model = model, psill = thetas[1], range = thetas[2], nugget = 1 - thetas[1])
  A <- variogramLine(vgmodel, dist_vector = D, covariance = TRUE)
  cholA <- try(chol(A), silent = TRUE)
  if (inherits(cholA, "try-error")) {
    return(Inf)
  } else {
    I <- FIM(A, D, model, thetas, perturbation) #Fisher information matrix
    cholI <- try(chol(I), silent = TRUE)
    if (inherits(cholI, "try-error")) {
      return(Inf)
      } else {
        invI <- solve(I)
        #compute distance matrix and correlation matrix of sampling points used for prediction
        D <- spDists(psample)
        A <- variogramLine(vgm(model = model, psill = thetas[1], range = thetas[2], nugget = 1 - thetas[1]),
                             dist_vector = D, covariance = TRUE)
        #extend correlation matrix A with a column and row with ones (ordinary kriging)
        nobs <- length(psample)
        B <- matrix(data = 1, nrow = nobs + 1, ncol = nobs + 1)
        B[1:nobs, 1:nobs] <- A
        B[nobs + 1, nobs + 1] <- 0

        #compute matrix with correlations between evaluation points and sampling points used for prediction
        D0 <- spDists(x = esample, y = psample)
        A0 <- variogramLine(vgm(model = model, psill = thetas[1], range = thetas[2], nugget = 1 - thetas[1]),
                                dist_vector = D0, covariance = TRUE)
        b <- cbind(A0, 1)
        #compute perturbed correlation matrices (pA and pA0)
        pA <- pA0 <- pB <- pb <- list()
        for (i in seq_len(length(thetas))) {
          thetas.pert <- thetas
          thetas.pert[i] <- (1 + perturbation) * thetas[i]
          pA[[i]] <- variogramLine(vgm(model = model, psill = thetas.pert[1], range = thetas.pert[2], nugget = 1 - thetas.pert[1]),
                                       dist_vector = D, covariance = TRUE)
          pA0[[i]] <- variogramLine(vgm(model = model, psill = thetas.pert[1], range = thetas.pert[2], nugget = 1 - thetas.pert[1]),
                                  dist_vector = D0, covariance = TRUE)
          pB[[i]] <- B
          pB[[i]][1:nobs, 1:nobs] <- pA[[i]]
          pb[[i]] <- cbind(pA0[[i]], 1)
        }
        #compute perturbed kriging variances (pvar)
        var <- numeric(length = length(esample)) #kriging variance
        pvar <- matrix(nrow = length(esample), ncol = length(thetas)) #matrix with perturbed kriging variances
        for (i in seq_len(length(esample))) {
          l <- solve(B, b[i,])
          var[i] <- 1 - l[1:nobs] %*% A0[i, ] - l[nobs + 1]
          for (j in seq_len(length(thetas))) {
            pl <- solve(pB[[j]], pb[[j]][i, ])
            pvar[i, j] <- 1 - pl[1:nobs] %*% pA0[[j]][i, ] - pl[nobs + 1]
          }
        }
        #approximate partial derivatives of kriging variance to correlogram parameters
        dvar <- list()
        for (i in seq_len(length(thetas))) {
          dvar[[i]] <- (pvar[, i] - var) / (thetas[i] * perturbation)
        }
        #compute variance of kriging variance for evaluation points.
        VKV <- numeric(length = length(var))
        for (i in seq_len(length(thetas))) {
          for (j in seq_len(length(thetas))) {
            VKVij <- invI[i, j] * dvar[[i]] * dvar[[j]]
            VKV <- VKV + VKVij
          }
        }
        MVKV <- mean(VKV)
        return(MVKV)
      }
  }
}


#' Mean Augmented Kriging Variance
#'
#' @param points \code{\link{data.frame}} with coordinates of fixed and free sampling points
#' @param esample \code{\link[sp]{SpatialPoints}} used for evaluation.
#' @param model \pkg{gstat} model type of a priori semivariogram model.
#' @param thetas parameters of semivariogram model.
#' @param perturbation proportion of a semivariogram parameter value added to
#'   that value.
#'
#' @importFrom sp "coordinates<-" spDists
#' @importFrom gstat variogramLine vgm
#'
#' @return mean augmented kriging variance.
#'
#' @export
MAKV <- function(points, esample, model, thetas, perturbation = 0.01)  {
  validate(esample = esample, model, thetas, perturbation)
  points <- as.data.frame(points)
  nobs <- nrow(points)
  coordinates(points) <- ~x + y
  D <- spDists(points)
  vgmodel <- vgm(model = model, psill = thetas[1], range = thetas[2], nugget = 1 - thetas[1])
  A <- variogramLine(vgmodel, dist_vector = D, covariance = TRUE)
  pA <- dA <- list()
  for (i in seq_len(length(thetas))) {
    thetas.pert <- thetas
    thetas.pert[i] <- (1 + perturbation) * thetas[i]
    vgmodel.pert <- vgm(model = model, psill = thetas.pert[1], range = thetas.pert[2], nugget = 1 - thetas.pert[1])
    pA[[i]] <- variogramLine(vgmodel.pert, dist_vector = D, covariance = TRUE)
    dA[[i]] <- (pA[[i]] - A) / (thetas[i] * perturbation)
  }

  cholA <- try(chol(A), silent = TRUE)
  if (inherits(cholA, "try-error")) {
    return(Inf)
  } else {
    I <- FIM(A, D, model, thetas, perturbation) #Fisher Information Matrix
    cholI <- try(chol(I), silent = TRUE)
    if (inherits(cholI, "try-error")) {
      return(Inf)
      } else {
        invI <- solve(I)

        nrowB <- nobs + 1
        B <- matrix(data = 1, nrow = nrowB, ncol = nrowB)
        B[1:nobs, 1:nobs] <- A
        B[nrowB, nrowB] <- 0

        #compute matrix with covariances between evaluation points and sampling points
        D0 <- spDists(x = esample, y = points)
        vgmodel <- vgm(model = model, psill = thetas[1], range = thetas[2], nugget = 1 - thetas[1])
        A0 <- variogramLine(vgmodel, dist_vector = D0, covariance = TRUE)
        #compute pB and pb by extending pA and pA0 with ones
        pB <- pA0 <- pb <- list()
        for (i in seq_len(length(thetas))) {
          pB[[i]] <- B
          pB[[i]][1:nobs, 1:nobs] <- pA[[i]]

          thetas.pert <- thetas
          thetas.pert[i] <- (1 + perturbation) * thetas[i]
          vgmodel.pert <- vgm(model = model, psill = thetas.pert[1], range = thetas.pert[2], nugget = 1 - thetas.pert[1])
          pA0[[i]] <- variogramLine(vgmodel.pert, dist_vector = D0, covariance = TRUE)
          pb[[i]] <- cbind(pA0[[i]], 1)
        }

        L <- matrix(nrow = length(esample), ncol = nobs) #matrix with kriging weights
        pL <- array(dim = c(length(esample), length(points), length(thetas))) #array with perturbed kriging weights
        var <- numeric(length = length(esample)) #kriging variance
        for (i in seq_len(length(esample))) {
          b <- c(A0[i, ], 1)
          l <- solve(B, b)
          L[i, ] <- l[1:nobs]
          var[i] <- 1 - l[1:nobs] %*% A0[i, ] - l[-(1:nobs)]
          for (j in seq_len(length(thetas))) {
            pl <- solve(pB[[j]], pb[[j]][i, ])
            pL[i, , j] <- pl[1:nobs]
          }
        }

        dL <- list()
        for (i in seq_len(length(thetas))) {
          dL[[i]] <- (pL[, , i] - L) / (thetas[i] * perturbation)
        }
        #tausq: expectation of additional variance due to uncertainty in ML estimates of variogram parameters, see Eq. 5 Lark and Marchant 2018
        tausq <- numeric(length = length(esample))
        tausqk <- 0
        for (k in seq_len(length(esample))) {
          for (i in seq_len(length(dL))) {
            for (j in seq_len(length(dL))) {
              tausqijk <- invI[i, j] * t(dL[[i]][k, ]) %*% A %*% dL[[j]][k, ]
              tausqk <- tausqk + tausqijk
            }
          }
          tausq[k] <- tausqk
          tausqk <- 0
        }
        AKV <- var + tausq
        MAKV <- mean(AKV)
        return(MAKV)
      }
    }
  }


#' Mean Estimation Adjusted Criterion
#'
#' @param points \code{\link{data.frame}} with coordinates of fixed and free sampling points
#' @param esample \code{\link[sp]{SpatialPoints}} used for evaluation.
#' @param model \pkg{gstat} model type of a priori semivariogram model.
#' @param thetas parameters of semivariogram model.
#' @param perturbation proportion of a semivariogram parameter value added to
#'   that value.
#'
#' @importFrom sp "coordinates<-" spDists
#' @importFrom gstat variogramLine vgm
#'
#' @return Mean estimation adjusted criterion.
#'
#' @export
MEAC <- function(points, esample, model, thetas, perturbation = 0.01)  {
  validate(esample = esample, model, thetas, perturbation)
  points <- as.data.frame(points)
  nobs <- nrow(points)
  coordinates(points) <- ~x + y
  D <- spDists(points)
  vgmodel <- vgm(
    model = model, psill = thetas[1], range = thetas[2],
    nugget = 1 - thetas[1])
  A <- variogramLine(vgmodel, dist_vector = D, covariance = TRUE)
  pA <- dA <- list()
  for (i in seq_len(length(thetas))) {
    thetas.pert <- thetas
    thetas.pert[i] <- (1 + perturbation) * thetas[i]
    vgmodel.pert <- vgm(
      model = model, psill = thetas.pert[1], range = thetas.pert[2],
      nugget = 1 - thetas.pert[1])
    pA[[i]] <- variogramLine(vgmodel.pert, dist_vector = D, covariance = TRUE)
    dA[[i]] <- (pA[[i]] - A) / (thetas[i] * perturbation)
  }

  cholA <- try(chol(A), silent = TRUE)
  if (inherits(cholA, "try-error")) {
    return(Inf)
    } else {
      I <- FIM(A, D, model, thetas, perturbation) #Fisher Information Matrix
      cholI <- try(chol(I), silent = TRUE)
      if (inherits(cholI, "try-error")) {
        return(Inf)
        } else {
          invI <- solve(I)

          nrowB <- nobs + 1
          B <- matrix(data = 1, nrow = nrowB, ncol = nrowB)
          B[1:nobs, 1:nobs] <- A
          B[nrowB, nrowB] <- 0

          #compute matrix with covariances between prediction nodes and sampling points
          D0 <- spDists(x = esample, y = points)
          vgmodel <- vgm(model = model, psill = thetas[1], range = thetas[2], nugget = 1 - thetas[1])
          A0 <- variogramLine(vgmodel, dist_vector = D0, covariance = TRUE)
          #compute pB and pb by extending pA and pA0 with ones
          pB <- pA0 <- pb <- list()
          for (i in seq_len(length(thetas))) {
            pB[[i]] <- B
            pB[[i]][1:nobs, 1:nobs] <- pA[[i]]

            thetas.pert <- thetas
            thetas.pert[i] <- (1 + perturbation) * thetas[i]
            vgmodel.pert <- vgm(model = model, psill = thetas.pert[1], range = thetas.pert[2], nugget = 1 - thetas.pert[1])
            pA0[[i]] <- variogramLine(vgmodel.pert, dist_vector = D0, covariance = TRUE)
            pb[[i]] <- cbind(pA0[[i]], 1)
          }

          L <- matrix(nrow = length(esample), ncol = nobs) #matrix with kriging weights
          pL <- array(dim = c(length(esample), length(points), length(thetas))) #array with perturbed kriging weights
          var <- numeric(length = length(esample)) #kriging variance
          pvar <- matrix(nrow = length(esample), ncol = length(thetas)) #matrix with perturbed kriging variances
          for (i in seq_len(length(esample))) {
            b <- c(A0[i, ], 1)
            l <- solve(B, b)
            L[i, ] <- l[1:nobs]
            var[i] <- 1 - l[1:nobs] %*% A0[i, ] - l[-(1:nobs)]
            for (j in seq_len(length(thetas))) {
              pl <- solve(pB[[j]], pb[[j]][i, ])
              pL[i, , j] <- pl[1:nobs]
              pvar[i, j] <- 1 - pl[1:nobs] %*% pA0[[j]][i, ] - pl[-(1:nobs)]
            }
          }

          dvar <- dL <- list()
          for (i in seq_len(length(thetas))) {
            dvar[[i]] <- (pvar[, i] - var) / (thetas[i] * perturbation)
            dL[[i]] <- (pL[, , i] - L) / (thetas[i] * perturbation)
          }
          #tausq: expectation of additional variance due to uncertainty in ML estimates of variogram parameters
          tausq <- numeric(length = length(esample))
          tausqk <- 0
          for (k in seq_len(length(esample))) {
            for (i in seq_len(length(dL))) {
              for (j in seq_len(length(dL))) {
                tausqijk <- invI[i, j] * t(dL[[i]][k, ]) %*% A %*% dL[[j]][k, ]
                tausqk <- tausqk + tausqijk
              }
            }
            tausq[k] <- tausqk
            tausqk <- 0
          }
          AKV <- var + tausq

          #VKV: variance of kriging variance
          VKV <- numeric(length = length(var))
          for (i in seq_len(length(dvar))) {
            for (j in seq_len(length(dvar))) {
              VKVij <- invI[i, j] * dvar[[i]] * dvar[[j]]
              VKV <- VKV + VKVij
            }
          }
          MEAC <- mean(AKV + VKV / (2 * var))
          MEAC
        }
    }
  }
