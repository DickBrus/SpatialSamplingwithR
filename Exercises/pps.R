#load  data
load(file="../data/Amazonia_5km.RData")

lnSWIR_neg <- gridAmazonia$lnSWIR2*-1
min_lnSWIR_neg <- min(lnSWIR_neg)
size <- lnSWIR_neg - min_lnSWIR_neg + 1

#check whether size variable is strictly positive (> 0)
summary(size)

#set sample size (number of draws)
n <- 100
p <- size/sum(size)

#Repeat pps sampling 10000 times
mz <- v_mz <- numeric(length=10000)
N <- nrow(gridAmazonia)
set.seed(314)

for (i in 1:10000) {
  mysample <- sample.int(N, size=n, replace=TRUE, prob=p)
  zexpand <- gridAmazonia$AGB[mysample]/p[mysample]
  tz <- mean(zexpand)
  mz[i] <- tz/N
  v_tz <- var(zexpand)/n
  v_mz[i] <- v_tz/N^2
}

hist(mz)
print(Vp_mz <- var(mz))
print(Ep_v_mz <- mean(v_mz))

#Compute true sampling variance of estimated mean with simple random sampling with replacement
print(v_mz_SI_true <-var(gridAmazonia$AGB)/d)

print(gain <- v_mz_SI_true/Vp_mz)
