library(spcosa)

s1 <- s2 <- 1:40 - 0.5
grid <- expand.grid(s1, s2)
names(grid) <- c("s1", "s2")

#compute distance of grid points to centre of square
d <- sqrt((grid$s1 - mean(grid$s1))^2 + (grid$s2 - mean(grid$s2))^2)

#select grid nodes with distance to centre smaller or equal to 20
plt <- grid[(d <= 20), ]

#compute compact geostrata
gridded(plt) <- ~ s1 + s2
set.seed(314)
myStrata <- spcosa::stratify(plt, nStrata = 6, equalArea = FALSE, nTry = 500)
#select the centres of the geostrata
mySample <- spcosa::spsample(myStrata)
#plot geostrata and sampling points (centres of geostrata)
plot(myStrata, mySample)

#save(myStrata,mySample, file = "../results/SpatialCoverageCircularPlot_6pnts.rda")

#check size of strata
myStrata <- as(myStrata, "data.frame")
table(myStrata$stratumId)
