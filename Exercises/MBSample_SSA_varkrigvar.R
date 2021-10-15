library(spsann)
library(sp)
library(matrixcalc)
library(ggplot2)


source ("../Rscripts/ObjectiveFunctions4MBSamplingVariogram.R")

load(file="../data/HunterValley.RData")

grd <- grdHunterValley
gridded(grd) <- ~Easting+Northing
candi <- spsample(grd, type="regular", cellsize=c(50,50))
candi <- as.data.frame(candi)
names(candi) <- c("x","y")

xi <- 0.5
phi <- 200
thetas <- c(xi,phi)

xgrid <- c(0,1,2,3)
ygrid <- xgrid
grid <- expand.grid(xgrid,ygrid)
names(grid) <- c("x","y")
spacing <- 300
grid$x <- grid$x*spacing
grid$y <- grid$y*spacing
coordinates(grid) <- ~x+y
myevalsample <- data.frame(x=450, y=450)
coordinates(myevalsample) <- ~x+y

schedule <- scheduleSPSANN(
  initial.acceptance=0.8,
  initial.temperature=0.0004,
  temperature.decrease=0.9,
  chains=300,
  chain.length=2,
  stopping=10,
  x.min=50,y.min=50,
  cellsize=50)

set.seed(314)

rslt <- optimUSER(
  points=100,
  candi=candi,
  fun=varkrigvar,
  grid=grid,
  esample=myevalsample,
  model="Exp",
  thetas=thetas,
  perturbation=0.01,
  schedule=schedule,
  track=TRUE)


mysample <- candi[rslt$points$id,]
#save(mysample, file="results/MBSample_varkrigvar_phi200nug05.RData")

ggplot(data=grdHunterValley) +
  geom_raster(mapping = aes(x=Easting/1000, y=Northing/1000),fill="grey")+
  geom_point(data = mysample, mapping=aes(x=x/1000, y=y/1000), shape=2, size=2)+
  scale_x_continuous(name = "Easting (km)")+
  scale_y_continuous(name = "Northing (km)")+
  coord_fixed() 
