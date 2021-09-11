library(sampling)
load("data/Voorst.RData")

#compute total number of pixels per stratum and stratum weights (relative size)
N_h <- tapply(grdVoorst$z, INDEX=grdVoorst$newstratum, FUN=length)
w_h <- N_h/sum(N_h)

#total sample size
n <- 40

#compute stratum sample sizes for proportional allocation 
n_h <- round(n*w_h)
sum(n_h)

#reduce sample size of largest stratum
print(n_h)
n_h[1] <- n_h[1] - 1

#select stratified simple random sample 
ord <- unique(grdVoorst$stratum)
set.seed(314)
units <- sampling::strata(grdVoorst, stratanames="stratum", size=n_h[ord], method="srswr")
mysample <- getdata(grdVoorst, units)

#estimate population means and standard error of estimator
mz_h <- tapply(mysample$z, INDEX=mysample$stratum, FUN=mean)
mz <- sum(w_h*mz_h)
S2z_h <- tapply(mysample$z, INDEX=mysample$stratum, FUN=var)
v_mz_h <- S2z_h/n_h
se_mz <- sqrt(sum(w_h^2*v_mz_h))



#compute true standard error
S2z_h_pop <- tapply(grdVoorst$z, INDEX=grdVoorst$stratum, FUN=var)
v_mz_h_true <- S2z_h_pop/n_h
print(se_mz_true <- sqrt(sum(w_h^2*v_mz_h_true)))



#compute 95% confidence interval estimate of population mean with function confint of package survey
library(survey)

#first add stratum weights to the selected sample
labels <- sort(unique(mysample$stratum))
lut <- data.frame(stratum=labels, weight=N_h/n_h)
mysample <- merge(x=mysample, y=lut)

#specify the sampling design
design_stsi <- svydesign(id=~1, strata=~stratum, weight=~weight, data=mysample)

#estimate the population mean and standard error of estimator
(res <- svymean(~z, design_stsi))

#approximate degrees of freedom
df_stsi <- degf(design_stsi)

#compute confidence interval estimate
confint(res, df=df_stsi, level=0.95)
