grdKandahar <- readRDS(file = "../data/grdKandahar.rds")

tz <- sum(grdKandahar$poppy)
tx <- sum(grdKandahar$agri)
(ratio <- tz / tx)

#fit homoscedastic model without intercept (residual variance constant)
model_homo <- lm(poppy ~ agri - 1, data = grdKandahar)
(slope_homo <- model_homo$coef)

#fit heteroscedastic model without intercept (residual variance proportional to covariate)
model_hetero <- lm(poppy ~ agri - 1, data = grdKandahar, weights = 1 / agri)
(slope_hetero <- model_hetero$coef)
