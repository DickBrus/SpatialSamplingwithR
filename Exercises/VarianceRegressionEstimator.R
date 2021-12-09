library(sswr)

grdAmazonia$lnSWIR2 <- log(grdAmazonia$SWIR2)
N <- nrow(grdAmazonia)

lm_exhaustive <- lm(AGB~lnSWIR2, data = grdAmazonia)
ab_pop <- coef(lm_exhaustive)
S2e_pop_exhaustive <- sum(lm_exhaustive$residuals^2) / nrow(grdAmazonia)

# set sample sizes
n <- c(10, 25, 100)

(Vp_mz_regr_exhaustive <- S2e_pop_exhaustive / n)

mx_pop <- mean(grdAmazonia$lnSWIR2)
Vp_mz_regr <- Ep_av_mz_regr <- numeric(length = 3)
mz_regr <- av_mz_regr <- numeric(length = 10000)

set.seed(314)
for (i in seq_len(length(n))) {
  for (j in 1:10000) {
    units <- sample(nrow(grdAmazonia), size = n[i], replace = FALSE)
    mysample <- grdAmazonia[units, ]
    mx_sam <- mean(mysample$lnSWIR2)
    mz_sam <- mean(mysample$AGB)

    #fit model on sample data
    lm_sam <- lm(AGB~lnSWIR2, data = mysample)
    ab_sam <- coef(lm_sam)

    #compute regression estimate
    mz_regr[j] <- mz_sam + ab_sam[2] * (mx_pop - mx_sam)
    S2e_sam <- sum(lm_sam$residuals^2) / (n[i] - 1)

    #approximate the variance of the regression estimator
    av_mz_regr[j] <- (1 - n[i] / N) * S2e_sam / n[i]
  }

Vp_mz_regr[i] <- var(mz_regr)
Ep_av_mz_regr[i] <- mean(av_mz_regr)
}

df <- data.frame(n, Vp_mz_regr_exhaustive, Vp_mz_regr, Ep_av_mz_regr)
#save(df, file = "results/VarRegressionEstimator.RData")
names(df) <- c("n", "VarMean_exh", "VarMean_sam", "ApproxVar")

#compute difference between variance of regression estimator with exhaustive fit of model and variance with sample fit of model,
#as a proportion of the variance with sample fit
with(df, (VarMean_exh - VarMean_sam) / VarMean_sam)
#same for difference between approximated variance and variance with sample fit
with(df, (ApproxVar - VarMean_sam) / VarMean_sam)

library(tidyverse)
df_lf <- df %>% pivot_longer(., cols = c("VarMean_exh", "VarMean_sam", "ApproxVar"))

library(ggplot2)
ggplot(df_lf) +
  geom_point(mapping = aes(x = n, y = value, shape = name), size = 2) +
  scale_shape_manual(values = c(1, 0, 2), name = "") +
  scale_x_continuous(name = "Sample size") +
  scale_y_continuous(name = "Variance")