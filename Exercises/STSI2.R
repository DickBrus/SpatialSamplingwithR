library(stratification)
library(forcats)

grdVoorst <- readRDS("../data/grdVoorst.rds")
grdVoorst$newstratum <- fct_collapse(grdVoorst$stratum,SA=c("EA","PA"))

#compute total number of pixels per stratum and stratum weights (relative size)
N_h <- tapply(grdVoorst$z, INDEX=grdVoorst$newstratum, FUN=length)
w_h <- N_h/sum(N_h)

#total sample size
n <- 40

#compute stratum sample sizes for proportional allocation 
n_h <- round(n*w_h)
sum(n_h)

#too many points, reduce sample size of merged stratum SA
n_h[2] <- n_h[2]-1

#the sampling variance can be computed without error from the simulated population
S2z_h_pop <- tapply(grdVoorst$z, INDEX=grdVoorst$newstratum,FUN=var) 
print(v_mz_4strata_true <- sum(w_h^2*S2z_h_pop/n_h))

#compare with true sampling variance obtained with original strata (same sample size, proportional allocation)

N_h <- tapply(grdVoorst$z,INDEX=grdVoorst$stratum,FUN=length)
w_h <- N_h/sum(N_h)

n_h <- round(n*w_h)
sum(n_h)
n_h[1] <- n_h[1]-1

S2z_h_pop <- tapply(grdVoorst$z, INDEX=grdVoorst$stratum,FUN=var) 
print(v_mz_5strata_true <- sum(w_h^2*S2z_h_pop/n_h))

#conclusion: sampling design with four strata only marginally less precise than design with five strata