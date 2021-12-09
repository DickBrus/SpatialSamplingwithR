#Load relevant packages

library(sp)
library(ggplot2)
library(sswr)

#Define function for simple random sampling of pairs of points with separation distance h

SIpairs <- function(h, n, area) {
  topo <- as(getGridTopology(area), "data.frame")
  cellsize <- topo$cellsize[1]
  xy <- coordinates(area)
  dxy <- numeric(length = 2)
  xypnts1 <- xypnts2 <- NULL
  i <- 1
  while (i <= n) {
    id1 <- sample(length(area), size = 1)
    xypnt1 <- xy[id1, ]
    xypnt1[1] <- jitter(xypnt1[1], amount = cellsize / 2)
    xypnt1[2] <- jitter(xypnt1[2], amount = cellsize / 2)
    angle <- runif(n = 1, min = 0, max = 2 * pi)
    dxy[1] <- h * sin(angle)
    dxy[2] <- h * cos(angle)
    xypnt2 <- xypnt1 + dxy
    xypnt2 <- as.data.frame(t(xypnt2))
    coordinates(xypnt2) <- ~s1 + s2
    inArea <- as.numeric(over(x = xypnt2, y = area))[1]
    if (!is.na(inArea)) {
      xypnts1 <- rbind(xypnts1, xypnt1)
      xypnts2 <- rbind(xypnts2, as.data.frame(xypnt2))
      i <- i + 1
      }
    rm(xypnt1, xypnt2)
  }
  cbind(xypnts1, xypnts2)
}

grd <- grdHunterValley
coordinates(grd) <- ~ s1 + s2
gridded(grd) <- TRUE

#Insert separation distances
h <- c(25, 50, 100, 200, 400)

#Select all point-pairs

set.seed(314)
samplesize <- 100

allpairs <- NULL
for (i in seq_len(length(h))) {
  pairs <- SIpairs(h = h[i], n = samplesize, area = grd)
  allpairs <- rbind(allpairs, pairs, make.row.names = FALSE)
}

#Overlay with grid

p1 <- allpairs[, c(1, 2)]
p2 <- allpairs[, c(3, 4)]
coordinates(p1) <- ~s1 + s2
z1 <- over(x = p1, y = grd)[4]
coordinates(p2) <- ~s1 + s2
z2 <- over(x = p2, y = grd)[4]

mysample <- data.frame(h = rep(h, each = samplesize), z1, z2)
names(mysample)[c(2, 3)] <- c("z1", "z2")
gammah <- vargammah <- numeric(length = length(h))

for (i in seq_len(length(h))) {
  units <- which(mysample$h == h[i])
  pairs.h <- mysample[units, ]
  gammah[i] <- mean((pairs.h$z1 - pairs.h$z2)^2, na.rm = TRUE) / 2
  vargammah[i] <- var((pairs.h$z1 - pairs.h$z2)^2, na.rm = TRUE) / (samplesize * 4)
}

#Plot sample semivariogram

samplevariogram <- data.frame(h, gammah, vargammah)
ggplot(data = samplevariogram) +
    geom_point(mapping = aes(x = h, y = gammah), size = 3) +
    scale_x_continuous(name = "Separation distance") +
    scale_y_continuous(name = "Semivariance", limits = c(0, NA))

#Fit model

sphericalnugget <- function(h, range, psill, nugget) {
    h <- h / range
    nugget + psill * ifelse(h < 1, (1.5 * h - 0.5 * h^3), 1)
}


fit.var <- nls(gammah~sphericalnugget(h, range, psill, nugget),
               data = samplevariogram,
               start = list(psill = 4, range = 200, nugget = 1),
               weights = 1 / vargammah,
               algorithm = "port",
               lower = c(0, 0, 0),
               trace = TRUE)

coef(fit.var)

#compute variance covariance matrix of estimated variogram parameters by bootstrapping

allpars <- NULL
nboot <- 100
for (j in 1:nboot) {
  #select bootstrap sample for each lag and compute semivariance
  gammah <- vargammah <- numeric(length = length(h))
  for (i in seq_len(length(h))) {
    units <- which(mysample$h == h[i])
    pairs <- mysample[units, ]
    mysampleunits <- sample(samplesize, size = samplesize, replace = TRUE)
    mybtpsample <- pairs[mysampleunits, ]
    gammah[i] <- mean((mybtpsample$z1 - mybtpsample$z2)^2, na.rm = TRUE) / 2
    vargammah[i] <- var((mybtpsample$z1 - mybtpsample$z2)^2, na.rm = TRUE) / (samplesize * 4)
  }

  #fit model
  samplevariogram <- data.frame(h, gammah, vargammah)
  tryCatch({
    fittedvariogram <- nls(gammah~sphericalnugget(h, range, psill, nugget),
                           data = samplevariogram,
                           start = list(psill = 4, range = 200, nugget = 1),
                           weights = 1 / vargammah,
                           algorithm = "port",
                           lower = c(0, 0, 0))
  pars <- coef(fittedvariogram)
  allpars <- rbind(allpars, pars)}, error = function(e){})
  }

#compute variance-covariance matrix of two variogram parameters
var(allpars)