library(spcosa)
library(sp)
library(matrixcalc)

# Annealing functions
source("../Rscripts/ObjectiveFunctions4MBSamplingVariogram.R")

gridded(grdHunterValley) <- ~ s1 + s2
candi <- spsample(grdHunterValley, type = "regular", cellsize = c(50, 50))
candi <- as.data.frame(candi)
names(candi) <- c("x", "y")

s <- 0.8 #ratio of spatial dependence c1/(c0+c1)
range <- 200
thetas <- c(s, range)

# Select spatial coverage sample for prediction. These locations are fixed, i.e. their locations are not optimised in simulated annealing
# Note that a spatial coverage sample is not strictly needed! The alternative is to optimise the coordinates of all points in SSA

# Choose number of units of spatial coverage sample
n  <- 80
set.seed(314)
myStrata <- stratify(grdHunterValley, nStrata = n, equalArea = FALSE, nTry = 10)
mySCsample <- as(spsample(myStrata), "SpatialPoints")

# Select initial supplemental sample
nsup <- 20
units <- sample(nrow(grdHunterValley), nsup)
mysupsample <- as(grdHunterValley[units, ], "SpatialPoints")

# Select evaluation sample
myevalsample <- spsample(x = grdHunterValley, n = 100, type = "regular", offset = c(0.5, 0.5))

# Set amount of perturbation of correlogram parameters
perturbation <- 0.01

schedule <- scheduleSPSANN(
  initial.acceptance = 0.8,
  initial.temperature = 0.008,
  temperature.decrease = 0.8,
  chains = 500,
  chain.length = 5,
  stopping = 10,
  x.min = 0, y.min = 0,
  cellsize = 50)

fixed <- coordinates(mySCsample)
names(fixed) <- c("x", "y")
pnts <- list(fixed = fixed, free = nsup)

set.seed(314)
res <- optimUSER(
  points = pnts,
  candi = candi,
  fun = MEAC,
  esample = myevalsample,
  model = "Exp",
  thetas = thetas,
  perturbation = 0.01,
  schedule = schedule,
  track = TRUE)

mysample <- res$points
saveRDS(res, file = "../results/MBSample_MEAC_phi200nug02_20sup_HunterValley.rds")
MEACopt <- tail(res$objective$energy$obj, 1)

units <- which(mysample$free == 1)
mysupsample <- mysample[units, c("x", "y")]

#compute shortest distance from supplemental point to spatial coverage points
mySCsampledf <- as(mySCsample, "data.frame")
D <- spDists(x = as(mysupsample, "matrix"), y = as(mySCsampledf, "matrix"))
DminAV <- apply(D, MARGIN = 1, FUN = min)

plot(myStrata) +
  geom_point(data = mySCsampledf, mapping = aes(x = s1, y = s2), shape = 1, size = 1.5) +
  geom_point(data = mysupsample, mapping = aes(x = x, y = y), shape = 2, size = 2, colour = "red")
