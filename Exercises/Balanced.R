library(sampling)
library(sswr)

# Select a 5 km x 5 km subgrid from grdAmazonia

gridded(grdAmazonia) <- ~ x1 + x2

subgrd <- spsample(
  x = grdAmazonia, type = "regular", cellsize = c(5000, 5000),
  offset = c(0.5, 0.5))

res <- over(subgrd, grdAmazonia)

subgrdAmazonia <- data.frame(coordinates(subgrd), res[,])

#compute natural logs of SWIR2
subgrdAmazonia$lnSWIR2 <- log(subgrdAmazonia$SWIR2)

#compute matrix with balancing variables (design matrix)
N <- nrow(subgrdAmazonia)
X <- cbind(rep(1, times = N), subgrdAmazonia$lnSWIR2, subgrdAmazonia$Terra_PP)

#compute equal inclusion probabilities for all N units
n <- 100
pi <- rep(n / N, times = N)

#select balanced sample
set.seed(314)
sample_ind <- samplecube(X = X, pik = pi, method = 1)
eps <- 1e-6
units <- which(sample_ind > (1 - eps))
mysample <- subgrdAmazonia[units, ]

#estimate the population mean of AGB by pi estimator
mz <- mean(mysample$AGB)

#function to estimate the population regression coefficients
estimate_b <- function(z, X, c) {
  cXX <- matrix(nrow = ncol(X), ncol = ncol(X), data = 0)
  cXz <- matrix(nrow = 1, ncol = ncol(X), data = 0)
  for (i in seq_len(length(z))) {
    x <- X[i, ]
    cXX_i <- c[i] * (x %*% t(x))
    cXX <- cXX + cXX_i
    cXz_i <- c[i] * t(x) * z[i]
    cXz <- cXz + cXz_i
  }
  b <- solve(cXX, t(cXz))
  return(b)
}

#estimate regression coefficients
pi <- rep(n / N, n)
c <- (1 - pi)
b <- estimate_b(z = mysample$AGB / pi, X = X[units, ] / pi, c = c)

#compute predicted AGB for all units
zpred <- X %*% b

#compute residuals for selected units
e <- mysample$AGB - zpred[units]

#estimate variance of pi estimator of population total
v_tz <- n / (n - ncol(X)) * sum(c * (e / pi)^2)

#finally estimate the standard error of the pi estimator of the population mean
se_mz <- sqrt(v_tz / N^2)

#Equation (9.7)
sqrt(sum(c * e^2) / (n * (n - 3)))
