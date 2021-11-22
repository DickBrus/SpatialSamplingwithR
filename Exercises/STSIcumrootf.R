library(sampling)
library(stratification)
library(ggplot2)

#load grid
load(file = "../results/Amazonia_5km.RData")

#set number of strata
H <- 10

#compute optimal strata
nclass <- H*100 #number of classes used in computing the histogram

#grdAmazonia must be sorted in ascending order by the columns used for optimal stratification, see help of function strata.cumrootf
grdAmazonia <- grdAmazonia[order(grdAmazonia$lnSWIR2),]

#NB n does not influence the optimal stratification
n <- 100
optstrata <- strata.cumrootf(x=grdAmazonia$lnSWIR2,n=n,Ls=H,nclass=nclass)

optstrata$bh #stratum boundaries
optstrata$Nh #stratum sizes

#add optimal strata to data frame
grdAmazonia$optstrata <- optstrata$stratumID

#compute stratum sample sizes for proportional allocation
w_h <- optstrata$Nh/sum(optstrata$Nh)
(n_h <- round(w_h*n))
sum(n_h)

#select stratified simple random sample
units <- sampling::strata(grdAmazonia,stratanames="optstrata",size=n_h,method="srswor")
mysample <- getdata(grdAmazonia,units)
mz_h <- tapply(mysample$AGB, INDEX=mysample$optstrata,FUN=mean)
(mz <- sum(w_h*mz_h))
S2z_h <- tapply(mysample$AGB, INDEX=mysample$optstrata,FUN=var)
(v_mz <- sum(w_h^2*(1-n_h/optstrata$Nh)*S2z_h/n_h))

#The sampling variance can be computed without error from the population
S2z_h_pop <- tapply(grdAmazonia$AGB, INDEX=grdAmazonia$optstrata,FUN=var) 
(v_mz_true <- sum(w_h^2*(1-n_h/optstrata$Nh)*S2z_h_pop/n_h))

#compute gain in precision due to stratification (stratification effect)
v_mz_SI_true <- (1-n/sum(optstrata$Nh))*var(grdAmazonia$AGB)/n
(gain <- v_mz_SI_true/v_mz_true)
