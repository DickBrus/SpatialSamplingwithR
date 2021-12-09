library(sswr)

#construct PSUs
w <- 500 #width of PSUs
s1bnd <- seq(min(grdVoorst$s1) + w, min(grdVoorst$s1) + 11  * w, w) + 12.5
s1f <- findInterval(grdVoorst$s1, s1bnd)
s2bnd <- min(grdVoorst$s2) + w + 12.5
s2f <- findInterval(grdVoorst$s2, s2bnd)
grdVoorst$psu <- as.character(interaction(s1f, s2f))

#Question 1

M <- nrow(grdVoorst)

n <- 4
m <- 10

mz_pop <- mean(grdVoorst$z)
tz_pop <- sum(grdVoorst$z)
mz_psu <- tapply(grdVoorst$z, INDEX = grdVoorst$psu, FUN = mean)
tz_psu <- tapply(grdVoorst$z, INDEX = grdVoorst$psu, FUN = sum)
M_psu <- tapply(grdVoorst$z, INDEX = grdVoorst$psu, FUN = length)
p <- M_psu  / M
S2z_psu <- tapply(grdVoorst$z, INDEX = grdVoorst$psu, FUN = var) #NB this is the variance of SSUs within the PSUs
S2z_within <- sum(p  * S2z_psu) #pooled within primary unit variance
S2mz_between <- sum(p  * (mz_psu - mz_pop)^2)

print(v_mz_true <- S2mz_between  / n + S2z_within  / (n * m))

##Note that this is equal to Eq. 11.33 in Cochran (1977), with z_i = M_psu  / M
print(v_mz_true <- 1  / M^2 *
        (sum(p  * (tz_psu  / p - tz_pop)^2)  / n +
           sum(M_psu^2  * S2z_psu  / p)  / (n  * m)))



#Question 3: optimisation of sample size given precision requirement

v.max <- 1 #maximum variance of estimated population mean

c1 <- 2
c2 <- 1

print(n <- 1  / v.max  * (sqrt(S2z_within  * S2mz_between)  * sqrt(c2  / c1) + S2mz_between))
print(m <- sqrt(S2z_within)  / sqrt(S2mz_between)  * sqrt(c1  / c2))



#Question 4: optimisation of sample size for a given budget

C.max <- 100
(n <- C.max  * sqrt(S2mz_between)  / (sqrt(S2z_within)  * sqrt(c1  * c2) + sqrt(S2mz_between)  * c1))




#Question 5: compute optimal sample sizes with package PracTools

library(PracTools)

#compute variance components
res <- BW2stagePPS(X = grdVoorst$z, pp = p, psuID = grdVoorst$psu)

#Compute optimal n and m given a precision requirement
#compute target coefficient of variation of estimated population total
CV0 <- sqrt(v.max  * M^2)  / tz_pop
clusOpt2(C1 = c1, C2 = c2, delta = res[6], unit.rv = res[3], k = res[5], CV0 = CV0, cal.sw = 2)

#Compute optimal n and m given a budget constraint
clusOpt2(C1 = c1, C2 = c2, delta = res[6], unit.rv = res[3], k = res[5], tot.cost = C.max, cal.sw = 1)

#Note that n and m are reversed:
#n as computed with PracTools is the number of selected secondary units per primary unit draw
#m is the number of primary unit draws
