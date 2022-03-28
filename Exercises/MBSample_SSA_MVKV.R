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

set.seed(314)
mystrata <- stratify(grdHunterValley, nStrata = 100, equalArea = FALSE, nTry = 10)
mysample_SC <- as(spsample(mystrata), "SpatialPoints")
mysample_eval <- spsample(x = grdHunterValley, n = 200, type = "regular", offset = c(0.5, 0.5))

xi <- 0.5
phi <- 200
thetas <- c(xi,phi)

schedule <- scheduleSPSANN(
  initial.acceptance = c(0.8, 0.95),
  initial.temperature = 0.002,
  temperature.decrease = 0.8,
  chains = 300,
  chain.length = 2,
  stopping = 10,
  x.min = 0, y.min = 0,
  cellsize = 50)

set.seed(314)
res <- optimUSER(
  points = 50,
  candi = candi,
  fun = MVKV,
  psample = mysample_SC,
  esample = mysample_eval,
  model = "Exp",
  thetas = thetas,
  perturbation = 0.01,
  schedule = schedule,
  track = TRUE)

mysample <- res$points
write_rds(mysample, file = "results/MBSample_MVKV_phi200nug05_HunterValley_50pnts.rds")

ggplot(data=as.data.frame(grdHunterValley)) +
  geom_raster(mapping = aes(x = s1 / 1000, y = s2 / 1000), fill="grey")+
  geom_point(data = mysample, mapping=aes(x = x / 1000, y = y / 1000), shape=2, size=2)+
  scale_x_continuous(name = "Easting (km)")+
  scale_y_continuous(name = "Northing (km)")+
  coord_fixed() 
