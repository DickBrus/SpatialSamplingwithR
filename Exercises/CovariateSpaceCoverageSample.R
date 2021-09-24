library(fields)
library(ggplot2)

# read data
load(file="../data/HunterValley.RData")

n <- 20
set.seed(314)
covars <- c("cti","ndvi")
myClusters <- kmeans(scale(grdHunterValley[,covars]), centers=n, iter.max=1000, nstart=40)
grdHunterValley$cluster <- myClusters$cluster

#Select locations closest to the centers of the clusters
res <- fields::rdist(x1=myClusters$centers, x2=scale(grdHunterValley[,covars]))
units <- apply(res, MARGIN=1, FUN=which.min)
myCSCsample <- grdHunterValley[units, covars]

ggplot(grdHunterValley) +
  geom_point(mapping=aes(x=ndvi,y=cti, colour=as.character(cluster)),alpha=0.5) +
  scale_colour_viridis_d() +
  geom_point(data=myCSCsample,mapping=aes(x=ndvi,y=cti), size=2.5, colour="red") +
  scale_x_continuous(name = "ndvi") +
  scale_y_continuous(name = "cti") +
  theme(legend.position="none")

#Question 2

set.seed(314)
covars <- c("cti","ndvi","elevation_m")
myClusters <- kmeans(scale(grdHunterValley[,covars]), centers=n, iter.max=1000, nstart=40)
grdHunterValley$cluster2 <- myClusters$cluster

#Select locations closest to the centers of the clusters
res <- fields::rdist(x1=myClusters$centers, x2=scale(grdHunterValley[,covars]))
units <- apply(res,MARGIN=1, FUN=which.min)
myCSCsample2 <- grdHunterValley[units, c("cti","ndvi")]

ggplot(grdHunterValley) +
  geom_point(mapping=aes(x=ndvi,y=cti, colour=as.character(cluster2)),alpha=0.5) +
  scale_colour_viridis_d() +
  geom_point(data=myCSCsample2, mapping=aes(x=ndvi,y=cti), size=2.5, colour="red") +
  scale_x_continuous(name = "ndvi") +
  scale_y_continuous(name = "cti") +
  theme(legend.position="none")

#save(myCSCsample,myCSCsample2, grdHunterValley, file="results/CSCsampling_HunterValley.RData")
