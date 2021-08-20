library(gstat)
library(rgdal)

vgmodel <- vgm(model="Sph", psill=966, range=44.6)

shpField <- readOGR(dsn="../data",layer="Leest5")
proj4string(shpField) <- NA_character_
mygrid <- spsample(shpField,type="regular", n=2000, offset=c(0.5,0.5)) %>% as(.,"data.frame")
H <- as.matrix(dist(mygrid))
G <- variogramLine(vgmodel, dist_vector=H)
m_semivar_field <- mean(G)
u <- qnorm(0.975)
lmax <- 20
ceiling(n_SI <- (u*sqrt(m_semivar_field)/(lmax/2))^2)

#required sample size for SY 
gridded(mygrid) <- ~x1+x2
n <- 5:40
Exi_SY <- numeric(length=length(n))
set.seed(314)
for (i in 1:length(n)) {
  m_semivar_SY <- numeric(length=10)
  for (j in 1:10) {
    mySYsample <- spsample(x=mygrid, n=n[i], type="regular") %>% as(.,"data.frame")
    H_SY <- as.matrix(dist(mySYsample))
    G_SY <- variogramLine(vgmodel, dist_vector=H_SY)
    m_semivar_SY[j] <- mean(G_SY)
  }
  Exi_SY[i] <- m_semivar_field - mean(m_semivar_SY)
}

#compute length of CI (alpha=0.05)
lengthCI <- u*sqrt(Exi_SY)*2

plot(x=n, y=lengthCI)

(df <- data.frame(n,lengthCI))
(n_SY <- min(n[which(lengthCI<lmax)]))

#design effect
deff <- Exi_SY/(m_semivar_field/n)
(df <- data.frame(n, deff))

