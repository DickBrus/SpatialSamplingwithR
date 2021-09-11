load("data/Amazonia_1km.RData")
gridAmazonia$lnSWIR2 <- log(gridAmazonia$SWIR2)

library(survey)
n <- 50
set.seed(314)
mx_pop <- mean(gridAmazonia$lnSWIR2) 

R <- 1000
mz <- se_mz <- mz_g <- se_mz_g <- numeric(length=R)
for (i in 1:R) {
  
  units <- sample.int(nrow(gridAmazonia), size=n, replace=FALSE)
  mysample <- gridAmazonia[units,c("AGB","lnSWIR2")]
  X <- matrix(nrow=n, ncol=2, data=1)
  X[,2] <- mysample$lnSWIR2
  XXinv <- solve(t(X) %*% X)
  Xz <- t(X) %*% mysample$AGB
  ab <- t(XXinv %*% Xz)
  lm_sample <- lm(AGB~lnSWIR2, data=mysample)

  mx_sam <- mean(mysample$lnSWIR2) 
  mz_sam <- mean(mysample$AGB)
  mz[i] <- mz_sam+ab[2]*(mx_pop-mx_sam)

  e <- residuals(lm_sample)
  S2e <- var(e)
  N <- nrow(gridAmazonia)
  se_mz[i] <- sqrt((1-n/N)*S2e/n)

  mysample$fpc <- N
  design_si <- svydesign(id=~1, data=mysample,  fpc=~fpc)
  populationtotals <- c(N,sum(gridAmazonia$lnSWIR2))
  mysample_cal <- calibrate(design_si, formula=~lnSWIR2, population=populationtotals, calfun="linear")

  out <- svymean(~AGB, mysample_cal)
  mz_g[i] <- coef(out)
  se_mz_g[i] <- SE(out)
}

all.equal(mz,mz_g)
df <- data.frame(se_mz,se_mz_g)
#save(df, file="results/SERegressionEstimator.RData")

summary(df)

