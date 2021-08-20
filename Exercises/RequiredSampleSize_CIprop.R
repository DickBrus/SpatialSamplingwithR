library(binomSamSize)
library(ggplot2)

#compute required sample sizes for a given p0 (prior proportion), as a function of d (half the width of CI)
p0 <- 0.1
d <- seq(from=0.05, to=0.45, by=0.05)
n_prop_wald <- numeric(length=length(d))
for (i in 1:length(d)){
  n_prop_wald[i] <-  ciss.wald(p0=p0, d=d[i], alpha=0.05)
}
plot(x=d, y=n_prop_wald)

#compute required sample sizes for a given d, as a function of p0
d <- 0.2
p0 <- seq(from=0.05, to=0.5, by=0.05)
n_prop_wald <- numeric(length=length(p0))
for (i in 1:length(p0)){
  n_prop_wald[i] <- ciss.wald(p0=p0[i], d=d, alpha=0.05)
}
plot(x=p0, y=n_prop_wald)

