library(gstat)
library(sf)

nugget <- 483
psill <- 966 - nugget #966 is the sill of the variogram without nugget
vgm_SphNug <- vgm(model = "Sph", nugget = nugget, psill = psill, range = 44.6)

#read shape (geopackage) file and make grid

field <- read_sf("../data", "Leest") %>%
  st_set_crs(NA_crs_)

mygrid <- st_make_grid(field, cellsize = 2, what = "centers")
mygrid <- mygrid[field] %>%
  st_coordinates(mygrid)

H <- as.matrix(dist(mygrid))
G <- variogramLine(vgm_SphNug, dist_vector = H)
#replace zeroes on diagonal by nugget
diag(G) <- nugget
m_semivar_field <- mean(G)
n <- 25
(Exi_VSI <- m_semivar_field / n)

#same result is obtained by computing matrix G with spherical variogram without nugget, and adding the nugget to the mean semivariance

vgm_Sph <- vgm(model = "Sph", psill = psill, range = 44.6)
G <- variogramLine(vgm_Sph, dist_vector = H)
m_semivar_field <- mean(G) + nugget
(Exi_VSI <- m_semivar_field / n)