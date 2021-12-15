library(spsann)
library(sp)
library(spcosa)
library(ggplot2)
library(sswr)

gridded(grdHunterValley) <- ~ s1 + s2
candi <- spsample(grdHunterValley, type="regular",
                  cellsize=c(50,50), offset = c(0.5, 0.5))
candi <- as.data.frame(candi)
names(candi) <- c("x","y")

n <- 100
set.seed(314)
myStrata <- stratify(grdHunterValley, nStrata=n, equalArea=FALSE, nTry=10)
mySCsample <- as(spsample(myStrata), "SpatialPoints")
myevalsample <- spsample(x=grdHunterValley, n=100, type="regular", offset=c(0.5,0.5))
grdHunterValley <- as(grdHunterValley, "data.frame")

xi <- 0.5
phi <- 200
thetas <- c(xi,phi)

schedule <- scheduleSPSANN(
  initial.acceptance=0.8,
  initial.temperature=0.0004,
  temperature.decrease=0.9,
  chains=300,
  chain.length=2,
  stopping=10,
  x.min=0, y.min=0,
  cellsize=50)


set.seed(314)
res <- optimUSER(
  points=100,
  candi=candi,
  fun = MVKV,
  psample=mySCsample,
  esample=myevalsample,
  model= "Exp",
  thetas=thetas,
  perturbation=0.01, 
  schedule=schedule,
  track=TRUE)

mysample <- candi[res$points$id,]
saveRDS(mysample, file = "results/MBSample_MVKV_phi200nug05_HunterValley.rds")

ggplot(data=grdHunterValley) +
  geom_raster(mapping = aes(x = s1 / 1000, y = s2 / 1000), fill="grey")+
  geom_point(data = mysample, mapping=aes(x = x / 1000, y = y / 1000), shape=2, size=2)+
  scale_x_continuous(name = "Easting (km)")+
  scale_y_continuous(name = "Northing (km)")+
  coord_fixed() 
