library(fields)
library(ggplot2)

# read data
grdHunterValley <- readRDS(file = "../data/grdHunterValley.rds")

n <- 20
set.seed(314)
covars <- c("cti", "ndvi", "elevation_m")
myClusters <- kmeans(scale(grdHunterValley[, covars]), centers = n, iter.max = 1000, nstart = 40)
grdHunterValley$cluster2 <- myClusters$cluster

#Select locations closest to the centres of the clusters
res <- fields::rdist(x1 = myClusters$centers, x2 = scale(grdHunterValley[, covars]))
units <- apply(res, MARGIN = 1, FUN = which.min)
myCSCsample <- grdHunterValley[units, ]

ggplot(grdHunterValley) +
  geom_point(mapping = aes(x = ndvi, y = cti, colour = as.character(cluster2)), alpha = 0.5) +
  scale_colour_viridis_d() +
  geom_point(data = myCSCsample, mapping = aes(x = ndvi, y = cti), size = 2.5, colour = "red") +
  scale_x_continuous(name = "ndvi") +
  scale_y_continuous(name = "cti") +
  theme(legend.position = "none")

#pdf(file = "GeometricSamples_HunterValley.pdf", width = 7, height = 3.5)
ggplot(data = grdHunterValley) +
  geom_raster(mapping = aes(x = s1 / 1000, y = s2 / 1000, fill = cti)) +
  geom_point(data = myCSCsample, mapping = aes(x = s1 / 1000, y = s2 / 1000), colour = "orange", size = 1) +
  scale_fill_continuous(name = "twi", type = "viridis") +
  scale_x_continuous(name = "Easting (km)") +
  scale_y_continuous(name = "Northing (km)") +
  coord_fixed()
#dev.off()