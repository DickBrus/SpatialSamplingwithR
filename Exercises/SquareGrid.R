
library(sp)

#convert to SpatialPixelsDataFrame
gridded(grdAmhara) <- ~s1 + s2

#select a square grid, either by specifying the number of grid points or the spacing.
#use argument offset so that the starting point of the grid is not selected randomly

#set required sample size (number of grid points)
n <- 100
mysample <- spsample(x = grdAmhara, n = n, type = "regular", offset = c(0.5, 0.5))

mysample <- as(mysample, "data.frame")

#compute sample size
(nrow(mysample))

#alternative: set spacing
spacing <- 10200
mysample2 <- spsample(
  x = grdAmhara, cellsize = spacing, type = "regular", offset = c(0.5, 0.5))
mysample2 <- as(mysample2, "data.frame")

#compute sample size
(nrow(mysample2))

#Write a for-loop to select 200 times a square grid with on average 100 points, with random start
set.seed(314)
samplesize <- numeric(length = 200)
for (i in 1:200) {
  mysample <- spsample(x = grdAmhara, n = n, type = "regular")
  samplesize[i] <- length(mysample)
}

summary(samplesize)
hist(samplesize)

#Now select a square grid of exactly 100 points.

for (i in 1:50) {
  mysample <- spsample(x = grdAmhara, n = n, type = "regular")
  if (length(mysample) == n) {
    break
    }
}
length(mysample)

#convert SpatialPoints to data.frame, and print the coordinates.

(mysample <- as(mysample, "data.frame"))
