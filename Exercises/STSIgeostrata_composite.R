# load packages
library(sp)
library(spcosa)
library(sswr)

# compute compact geographical strata
# set random seed (for reproduction of results)
set.seed(314)
#convert grdVoorst to SpatialPixelsDataFrame
gridded(grdVoorst) <- ~ s1 + s2
myStrata <- stratify(grdVoorst, nStrata = 20, equalArea = TRUE, nTry = 1)
plot(myStrata)

# select randomly n locations from the strata
n_h <- 4
#mysample <- spsample(myStrata, n = n_h, type = "composite")
mysample <- spcosa::spsample(myStrata, n = n_h)
plot(myStrata, mysample)

#overlay selected points with simulated field to extract study variable at sampling points
mysample <- as(mysample, "SpatialPoints")
mysampledata <- over(mysample, grdVoorst)
mysample_df <- data.frame(as(mysample, "data.frame"), z = mysampledata$z)

head(mysample_df)

#add column with identifier for composite
mysample_df$composite <- rep(1:n_h)

#compute means for composites
mz_composites <- tapply(mysample_df$z, INDEX = mysample_df$composite, FUN = mean)

#Estimate mean for study area
print(mz <- mean(mz_composites))

#Estimate sampling variance of estimated mean
print(v_mz <- var(mz_composites) / n_h)
