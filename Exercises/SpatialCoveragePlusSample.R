library(spcosa)
library(sp)

#Read data with coordinates and other attributes of fine grid (discretisation of study area)

load(file="../data/HunterValley.RData")
grd <- grdHunterValley

head(grd)
N <- nrow(grd)

# Choose number of locations of spatial coverage sample

n<-90

#Compute clusters (geostrata) and select centers

coordinates(grd) <- ~Easting+Northing
gridded(grd) <- TRUE
set.seed(314)
myStrata <- stratify(grd, nStrata = n, equalArea=FALSE, nTry=10)
mySCsample <- spsample(myStrata)

# Compute average distance between neighbouring points of spatial coverage sample
gridTopology <- as(getGridTopology(grd), "data.frame")
A <- N *  gridTopology$cellsize[1]^2
d <- sqrt(A/n)

# Specify separation distances and subsample sizes
h <- c(20)
m <- c(10)

# Select random subsample from the spatial coverage sample

mySCsample.df <- as(mySCsample, "data.frame")
set.seed(314)
units <- sample(nrow(mySCsample.df),sum(m))
mySubsample <- mySCsample.df[units,]

# Select locations in random direction at distances h from subsample

first <- 1
final <- cumsum(m)
plus <- NULL
for (i in 1:length(h)) {
  dxy <- matrix(nrow=m[i],ncol=2)
  angle <- runif(n=m[i],min=0,max=2*pi)
  dxy[,1] <- h[i]*sin(angle)
  dxy[,2] <- h[i]*cos(angle)
  plus.h <- mySubsample[c(first:final[i]),]+dxy
  plus <- rbind(plus,plus.h)
  first <- final[i]+1
}

# Plot strata, spatial coverage sample and supplemental sample

library(ggplot2)
pdf(file = "SpatialCoveragePlusSample.pdf", width = 7, height = 7)
plot(myStrata,mySCsample) +
    geom_point(data = plus, mapping = aes(x= Easting, y =Northing), shape =2,size=1.5 )
dev.off()
