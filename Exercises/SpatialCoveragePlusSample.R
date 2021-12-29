library(spcosa)
library(sp)
library(sswr)

grdHunterValley

N <- nrow(grdHunterValley)

# Choose number of locations of spatial coverage sample

n <- 90

#Compute clusters (geostrata) and select centres

gridded(grdHunterValley) <- ~ s1 + s2
set.seed(314)
myStrata <- stratify(grdHunterValley, nStrata = n, equalArea = FALSE, nTry = 10)
mySCsample <- spsample(myStrata)

# Compute average distance between neighbouring points of spatial coverage sample
grdHunterValley <- as_tibble(grdHunterValley)
cellsize <- min(diff(sort(unique(grdHunterValley$s1))))
A <- N *  cellsize^2
d <- sqrt(A / n)

# Specify separation distances and subsample sizes
h <- c(20)
m <- c(10)

# Select random subsample from the spatial coverage sample

mySCsample.df <- as(mySCsample, "data.frame")
set.seed(314)
units <- sample(nrow(mySCsample.df), sum(m))
mySubsample <- mySCsample.df[units, ]

# Select locations in random direction at distances h from subsample

first <- 1
final <- cumsum(m)
plus <- NULL
for (i in seq_len(length(h))) {
  dxy <- matrix(nrow = m[i], ncol = 2)
  angle <- runif(n = m[i], min = 0, max = 2 * pi)
  dxy[, 1] <- h[i] * sin(angle)
  dxy[, 2] <- h[i] * cos(angle)
  plus_h <- mySubsample[c(first:final[i]), ] + dxy
  plus <- rbind(plus, plus_h)
  first <- final[i] + 1
}

# Plot strata, spatial coverage sample and supplemental sample

library(ggplot2)
#pdf(file = "SpatialCoveragePlusSample.pdf", width = 7, height = 7)
plot(myStrata, mySCsample) +
    geom_point(data = plus,
               mapping = aes(x = s1, y = s2), shape = 2, size = 1.5)
#dev.off()
