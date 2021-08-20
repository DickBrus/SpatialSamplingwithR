load(file="../data/Voorst.RData")

n <- 40
N <- nrow(grdVoorst)
set.seed(314)
units <- sample.int(N, size = n, replace = TRUE)
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
print(ind <- (mz_pop >lower95 & mz_pop < upper95))

#compute 90% confidence interval
margin90 <- qt(0.95, n-1, lower.tail=TRUE)*se_mz
print(lower90 <- mz - margin90)
print(upper90 <- mz + margin90)

#Estimate total mass of soil organic matter in topsoil (0 - 20 cm) of Voorst

#Compute the mass (kg) of 1 ha of soil until a depth of 0.2 m;
#Use a bulk density of soil of 1.4 g/cm3
soilmass_1ha <-  1.4 * 0.2 * 10^7

#Compute total area in ha of Voorst; 1 pixel of Voorst is 0.25 * 0.25 ha
ha_pop <- N*0.25^2

#Compute total mass of SOM in Voorst; unit of z is g/kg soil
soilmass <- soilmass_1ha * ha_pop

#Estimate total mass of SOM
print(tz <- mz * soilmass)

#Estimate the standard error of estimated total
v_tz <- v_mz*soilmass^2
print(se_tz <- sqrt(v_tz))