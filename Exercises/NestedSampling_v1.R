# First version of nested sampling in which in each stage one point is selected
# at some chosen distance h from all points selected in previous stages

# Only for balanced nested samples!

library(sp)
library(gstat)
library(ggplot2)
library(tidyverse)

#Define function for random selection of point at some chosen distance from a starting point

SelectPoint <- function(start, h, area) {
  dxy <- numeric(length = 2)
  inArea <- NA
  while (is.na(inArea)) {
    angle <- runif(n = 1, min = 0, max = 2 * pi)
    dxy[1] <- h * sin(angle)
    dxy[2] <- h * cos(angle)
    xypnt <- start + dxy
    coordinates(xypnt) <- ~s1 + s2
    inArea <- as.numeric(over(x = xypnt, y = area))[1]
  }
  xypoint <- as.data.frame(xypnt)
  xypoint
}

#Make a copy of grdHunterValley. The copy is used later to select a main station.
#The data.frame grdHunterValley is converted to a SpatialPixelsDataFrame

grid <- grdHunterValley
gridded(grdHunterValley) <- c("s1", "s2")

#Define separation distances
lags <- c(1000, 500, 200, 100, 50)

#Select main station
set.seed(123)
id <- sample(nrow(grid), 1)
mainstation <- grid[id, c(1, 2)]

#Select randomly one  point at distance lag[1] from main station
newpnt <- SelectPoint(start = mainstation, h = lags[1], area = grdHunterValley)
allstarts <- rbind(mainstation, newpnt)

for (j in 2:length(lags)) {
  newpnts <- NULL
  for (i in seq_len(nrow(allstarts))) {
    pnts <- SelectPoint(start = allstarts[i, ], h = lags[j], area = grdHunterValley)
    newpnts <- rbind(newpnts, pnts)
  }
  allstarts <- rbind(allstarts, newpnts)
}

nestedsample <- allstarts

#Add factors to the data frame nestedsample

nestedsample$f1 <- rep(1:2)
nestedsample$f2 <- rep(1:2, each = 2)
nestedsample$f3 <- rep(1:2, each = 4)
nestedsample$f4 <- rep(1:2, each = 8)

#Overlay with grdHunterValley and add the fourth attribute (cti) to the data frame nestedsample

coordinates(nestedsample) <- ~s1 + s2
res <- (over(nestedsample, grdHunterValley))[4]
nestedsample <- as(nestedsample, "data.frame")
nestedsample$cti <- res$cti

#Fit hierarchical ANOVA model by REML and extract the variance components

library(nlme)
lmodel <- lme(cti ~ 1, data = nestedsample, random = ~ 1 | f1 / f2 / f3 / f4)
res <- as.matrix(VarCorr(lmodel))
sigmas <- rev(as.numeric(res[c(2, 4, 6, 8, 9), 1]))
(semivar <- cumsum(sigmas))

#Plot sample semivariogram

df <- data.frame(rev(lags), semivar)
names(df)[1] <- "distance"
ggplot(df) +
  geom_point(mapping = aes(x = distance, y = semivar)) +
  scale_x_continuous(name = "Distance (m)") +
  scale_y_continuous(name = "Semivariance", limits = c(0, NA))
