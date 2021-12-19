library(sswr)

grdAmazonia$lnSWIR2 <- log(grdAmazonia$SWIR2)

library(survey)
n <- 50
set.seed(314)
mx_pop <- mean(grdAmazonia$lnSWIR2)

R <- 1000
mz <- se_mz <- mz_g <- se_mz_g <- numeric(length = R)
for (i in 1:R) {

  units <- sample(nrow(grdAmazonia), size = n, replace = FALSE)
  mysample <- grdAmazonia[units, c("AGB", "lnSWIR2")]
  lm_sample <- lm(AGB~lnSWIR2, data = mysample)
#inclusion probabilities are equal, so we can estimate regression coefficients by OLS
  ab <- coef(lm_sample)

  mx_sam <- mean(mysample$lnSWIR2)
  mz_sam <- mean(mysample$AGB)
  mz[i] <- mz_sam + ab[2] * (mx_pop - mx_sam)

  e <- residuals(lm_sample)
  S2e <- var(e)
  N <- nrow(grdAmazonia)
  se_mz[i] <- sqrt((1 - n / N) * S2e / n)

  mysample$fpc <- N
  design_si <- svydesign(id = ~ 1, data = mysample,  fpc = ~ fpc)
  populationtotals <- c(N, sum(grdAmazonia$lnSWIR2))
  mysample_cal <- calibrate(design_si, formula = ~ lnSWIR2,
                            population = populationtotals, calfun = "linear")

  out <- svymean(~AGB, mysample_cal)
  mz_g[i] <- coef(out)
  se_mz_g[i] <- SE(out)
}

all.equal(mz, mz_g)
df <- data.frame(se_mz, se_mz_g)
#save(df, file = "results/SERegressionEstimator.rda")

summary(df)