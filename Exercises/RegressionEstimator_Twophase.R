grdAmazonia <- readRDS(file = "../data/grdAmazonia.rds")
grdAmazonia$lnSWIR2 <- log(grdAmazonia$SWIR2)

N <- nrow(grdAmazonia)
#set sample sizes
n1 <- 250
n2 <- 100

mz_reg2ph <- av_mz_reg2ph <- mz_regr <- numeric(length = 10000)

set.seed(314)

for (i in 1:10000) {
  #select large simple random sample
  units <- sample(nrow(grdAmazonia), size = n1, replace = FALSE)
  mysample <- grdAmazonia[units, ]

  #fit simple linear model
  lm_sam <- lm(AGB ~ lnSWIR2, data = mysample)

  #estimate the population mean of AGB with the regression estimator
  mz_regr[i] <- mean(mysample$AGB) +
    coef(lm_sam)[2] * (mean(grdAmazonia$lnSWIR2) - mean(mysample$lnSWIR2))

  #subsample the selected sample
  units_subsam <- sample(nrow(mysample), size = n2, replace = FALSE)
  mysubsample <- mysample[units_subsam, ]

  #fit simple linear model on subsample
  lm_subsam <- lm(AGB~lnSWIR2, data = mysubsample)

  #compute regression estimator for two-phase sampling
  mz_reg2ph[i] <- mean(mysubsample$AGB) +
    coef(lm_subsam)[2] * (mean(mysample$lnSWIR2) - mean(mysubsample$lnSWIR2))

  #approximate the sampling variance of the regression estimator
  av_mz_reg2ph[i] <- (1 - n1 / N) * var(mysubsample$AGB) / n1 +
    (1 - n2 / n1) * (sum(lm_subsam$residuals^2) / (n2 - 2))/n1
}

(var(mz_regr))
(var(mz_reg2ph))
(mean(av_mz_reg2ph))

#save(mz_regr, mz_reg2ph, av_mz_reg2ph, file = "results/RegressionEstimator_Twophase.RData")
