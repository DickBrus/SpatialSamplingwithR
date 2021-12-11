library(tidyverse)
library(sf)
library(spcosa)

# set random seed (for reproduction of results)
set.seed(31415)

# read field of interest 
shpField <- read_sf("dat", "Leest") %>%
  st_set_crs(NA_crs_)

#shpField <- read_sf("dat", "Leest") %>%
#  st_set_crs(32631)

# compute compact geographical strata; either use argument cellSize or nGridCells
myStrata <- shpField %>%
  as_Spatial %>%
  stratify(nStrata = 10, cellSize = 2, equalArea = TRUE, nTry = 3)

plot(myStrata)

# obtain the areas of the strata
print(A_h <- getArea(myStrata))

# select randomly n points from the strata
mySample <- spsample(myStrata, n = 2)
plot(myStrata, mySample)

# convert mySample to SpatialPoints
samplingPoints <- as(mySample, "SpatialPoints")


# of naar sf
# samplingPoints <- mySample %>%
#     as("SpatialPoints") %>%
#     st_as_sf
