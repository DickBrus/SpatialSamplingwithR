load(file="../data/Voorst.RData")
load(file="data/Voorst.RData")

n <- 40
N <- nrow(grdVoorst)
set.seed(314)
units <- sample(N, size = n, replace = TRUE)
mysample <- grdVoorst[units,]
print(mz <- mean(mysample$z))

v_mz <- var(mysample$z)/n
print(se_mz <- sqrt(v_mz))

#compute lower and upper bound of 95% confidence interval using the t distribution
margin95 <- qt(0.975, n-1, lower.tail=TRUE)*se_mz
print(lower95 <- mz - margin95)
print(upper95 <- mz + margin95)

#check whether 95% confidence interval covers the true mean
print(mz_pop <- mean(grdVoorst$z))
print(ind <- (mz_pop > lower95 & mz_pop < upper95))

#compute 90% confidence interval
margin90 <- qt(0.95, n-1, lower.tail=TRUE)*se_mz
print(lower90 <- mz - margin90)
print(upper90 <- mz + margin90)

#Estimate total mass of soil organic matter in topsoil (0 - 30 cm) of Voorst

#Compute total volume of soil in layer 0 - 30 cm, in m3. 1 pixel of Voorst is 25 m * 25 m
vol_soil <- N * 25^2  * 0.3

#Compute total mass of soil in kg. Use bulk density of 1.5 g/cm3 = 1500 kg/m3
mass_soil <- vol_soil * 1500

#Estimate total mass of SOM in Mg. Note the units of mz is g/kg
print(tz <- mz * mass_soil * 10^-6)

#Estimate the standard error of estimated total
v_tz <- v_mz * mass_soil^2 *10^-12
print(se_tz <- sqrt(v_tz))
