library(sp)
library(spcosa)
library(ggplot2)
library(terra)

#read data frame with coordinates (and other attributes) of fine grid (discretisation of study area)
rmap <- rast(x="data/Elevation_Xuancheng.tif")
rmap
rmap <- as.data.frame(rmap, xy = TRUE, na.rm = TRUE)

gridded(rmap) <- ~ x + y
subgrid <- spsample(rmap, type = "regular", cellsize = c(900, 900), offset = c(0.5, 0.5))
length(subgrid)
subgrid <- as(subgrid, "data.frame")
rmap <- as(rmap, "data.frame")

#load existing sampling points
legacy <- sampleXuancheng[sampleXuancheng$sample == "iPSM", ]

#plot prior points
ggplot(rmap) +
  geom_raster(mapping = aes(x = x / 1000, y = y / 1000), fill = "grey") +
  geom_point(data = legacy, mapping = aes(x = s1 / 1000, y = s2 / 1000), size = 2) +
  scale_x_continuous(name = "Easting (km)") +
  scale_y_continuous(name = "Northing (km)") +
  coord_fixed(ratio = 1)

#convert data.frame subgrid to SpatialPixelsDataFrame
gridded(subgrid) <- ~ x1 + x2

#set number of new sampling points to be selected
n <- 100

#compute total sample size (existing points + new points)
ntot <- n + nrow(legacy)

#convert data.frame legacy (existing points) to SpatialPoints
legacy <- SpatialPoints(coords = cbind(legacy$s1, legacy$s2))

#compute geostrata with argument priorPoints
set.seed(314)
mystrata <- spcosa::stratify(subgrid, nStrata = ntot, priorPoints = legacy, nTry = 10)

#select sampling points of infill sample (centers of geostrata)
mysample <- spcosa::spsample(mystrata)

#plot geostrata and sampling points (centers of geostrata)
plot(mystrata, mysample)

#select the new points from mysample
units <- which(mysample@isPriorPoint == FALSE)

#convert mysample to data.frame
mysample <- as(mysample, "data.frame")
mysamplenew <- mysample[units, ]
