library(spsann)
library(sp)
library(matrixcalc)
library(ggplot2)

source("../Rscripts/ObjectiveFunctions4MBSamplingVariogram.R")

gridded(grdHunterValley) <- ~ s1 + s2
candi <- spsample(grdHunterValley, type = "regular", cellsize = c(50, 50))
candi <- as.data.frame(candi)
names(candi) <- c("x", "y")

xi <- 0.5
phi <- 200
thetas <- c(xi, phi)

schedule <- scheduleSPSANN(
  initial.acceptance = 0.8,
  initial.temperature = 0.04,
  temperature.decrease = 0.9,
  chains = 300,
  chain.length = 2,
  stopping = 10,
  x.min = 0, y.min = 0,
  cellsize = 50)

set.seed(314)
rslt <- optimUSER(
  points = 100,
  candi = candi,
  fun = logdet,
  model = "Exp",
  thetas = thetas,
  perturbation = 0.01,
  schedule = schedule,
  track = TRUE)

mysample <- candi[rslt$points$id, ]

saveRDS(mysample, file = "results/MBSample_logdet_phi200nug05_HunterValley.rds")

ggplot(data = grdHunterValley) +
  geom_raster(mapping = aes(x = s1 / 1000, y = s2 / 1000), fill = "grey") +
  geom_point(data = mysample, mapping = aes(x = x / 1000, y = y / 1000), shape = 2, size = 2) +
  scale_x_continuous(name = "Easting (km)") +
  scale_y_continuous(name = "Northing (km)") +
  coord_fixed()
