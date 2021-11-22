library(gstat)
nugget <- 483
psill <- 966 - nugget #966 is the sill of the variogram without nugget
vgm_SphNug <- vgm(model="Sph", nugget=nugget, psill=psill, range=44.6)

library(rgdal)
shpField <- readOGR(dsn="../data",layer="Leest")
proj4string(shpField) <- NA_character_
mygrid <- spsample(shpField,type="regular", n=2000, offset=c(0.5,0.5)) %>% as(.,"data.frame")
H <- as.matrix(dist(mygrid))
G <- variogramLine(vgm_SphNug, dist_vector=H)
#replace zeroes on diagonal by nugget
diag(G) <- nugget
m_semivar_field <- mean(G)
n <- 25
(Exi_VSI <- m_semivar_field/n)

#same result is obtained by computing matrix G with spherical variogram without nugget, and adding the nugget to the mean semivariance

vgm_Sph <- vgm(model="Sph", psill=psill, range=44.6)
G <- variogramLine(vgm_Sph,dist_vector=H)
m_semivar_field <- mean(G)+nugget
(Exi_VSI <- m_semivar_field/n)




