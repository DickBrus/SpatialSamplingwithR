library(gstat)
library(sf)

#specify semivariogram model as gstat object
vgmodel <- vgm(model = "Sph", psill = 966, range = 44.6)

#read geopackage file and make grid

field <- read_sf("../data/Leest.gpkg") %>%
  st_set_crs(NA_crs_)

mygrid <- st_make_grid(field, cellsize = 2, what = "centers")
mygrid <- mygrid[field] %>%
  st_coordinates(mygrid)



#compute required sample size for SI


H <- as.matrix(dist(mygrid))
G <- variogramLine(vgmodel, dist_vector = H)
m_semivar_field <- mean(G)
u <- qnorm(0.975)
lmax <- 20
ceiling(n_SI <- (u * sqrt(m_semivar_field) / (lmax / 2))^2)



#compute required sample size for SY


#convert mygrid to tibble (data.frame) and next to SpatialPixels object
mygrid <- mygrid %>%
  as_tibble() %>%
  setNames(c("x1", "x2"))
gridded(mygrid) <- ~ x1 + x2

n <- 5:40
Exi_SY <- numeric(length = length(n))
set.seed(314)
for (i in seq_len(length(n))) {
  m_semivar_SY <- numeric(length = 10)
  for (j in 1:10) {
    mySYsample <- spsample(x = mygrid, n = n[i], type = "regular") %>% 
      as("data.frame")
    H_SY <- as.matrix(dist(mySYsample))
    G_SY <- variogramLine(vgmodel, dist_vector = H_SY)
    m_semivar_SY[j] <- mean(G_SY)
  }
  Exi_SY[i] <- m_semivar_field - mean(m_semivar_SY)
}

#compute length of CI (alpha = 0.05)
lengthCI <- u * sqrt(Exi_SY) * 2

plot(x = n, y = lengthCI)

(df <- data.frame(n, lengthCI))
(n_SY <- min(n[which(lengthCI < lmax)]))

#design effect
deff <- Exi_SY / (m_semivar_field / n)
(df <- data.frame(n, deff))
