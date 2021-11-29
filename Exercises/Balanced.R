library(sampling)

grdAmazonia <- readRDS(file = "../results/grdAmazonia_5km.rds")

#compute natural logs of SWIR2
grdAmazonia$lnSWIR2 <- log(grdAmazonia$SWIR2)

#compute matrix with balancing variables (design matrix)
N <- nrow(grdAmazonia)
X <- cbind(rep(1, times = N), grdAmazonia$lnSWIR2, grdAmazonia$Terra_PP)

#compute equal inclusion probabilities for all N units
n <- 100
pi <- rep(n / N, times = N)

#select balanced sample
set.seed(314)
sample_ind <- samplecube(X = X, pik = pi, method = 1)
eps <- 1e-6
units <- which(sample_ind > (1 - eps))
mysample <- grdAmazonia[units, ]

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
