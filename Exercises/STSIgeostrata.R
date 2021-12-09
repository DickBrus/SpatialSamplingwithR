# load packages
library(spcosa)
library(rgdal)
library(sswr)

# set random seed (for reproduction of results)
set.seed(31415)

# read field of interest
shpField <- readOGR(dsn = "data", layer = "Leest", verbose = FALSE)

#shpField <- readOGR(dsn = "../data", layer = "Leest", verbose = FALSE)
proj4string(shpField) <- NA_character_

# compute compact geographical strata; either use argument cellSize or nGridCells
myStrata <- spcosa::stratify(
  shpField, nStrata = 10, cellSize = 2, equalArea = TRUE, nTry = 3)
#myStrata <- stratify(shpField, nStrata = 10, nGridCells = 2500, equalArea = TRUE, nTry = 3)
plot(myStrata)

# obtain the areas of the strata
print(A_h <- getArea(myStrata))

# select randomly n points from the strata
mySample <- spsample(myStrata, n = 2)
plot(myStrata, mySample)

# convert mySample to SpatialPoints
samplingPoints <- as(mySample, "SpatialPoints")
