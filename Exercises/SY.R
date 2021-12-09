library(sp)
library(ggplot2)
library(sswr)

#convert to SpatialPixelsDataFrame
gridded(grdVoorst) <- ~s1 + s2

#specify expected sample size
n <- 20
set.seed(314)

#select a first systematic random sample
mysamplecrds <- spsample(x = grdVoorst, n = n, type = "regular")
#overlay selected points with grdVoorst to extract values of z
mysample <- over(mysamplecrds, grdVoorst)
mysample_1 <- data.frame(coordinates(mysamplecrds), mysample)

#select a second systematic random sample
mysamplecrds <- spsample(x = grdVoorst, n = n, type = "regular")
#overlay selected points with grdVoorst to extract values of z
mysample <- over(mysamplecrds, grdVoorst)
mysample_2 <- data.frame(coordinates(mysamplecrds), mysample)

ggplot() +
  geom_raster(data = as(grdVoorst, "data.frame"),
              mapping = aes(x = s1 / 1000, y = s2 / 1000), fill = "grey") +
  geom_point(data = mysample_1,
             mapping = aes(x = x1 / 1000, y = x2 / 1000), size = 2) +
  geom_point(data = mysample_2,
             mapping = aes(x = x1 / 1000, y = x2 / 1000), size = 2, colour = "red") +
  scale_x_continuous(name = "Easting (km)") +
  scale_y_continuous(name = "Northing (km)") +
  coord_fixed()

mz <- numeric(length = 2)
mz[1] <- mean(mysample_1$z)
mz[2] <- mean(mysample_2$z)
print(mz)

print(v_mz <- var(mz))

#compute ultimate estimate of population mean and its sampling variance
print(mz_final <- mean(mz))

print(v_mz_final <- v_mz / 4)
