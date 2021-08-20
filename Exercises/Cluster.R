load(file="../data/Voorst.RData")

#construct zones
cellsize <- 25 #size of pixels
w <- 1000 #width of zones
s1bnd <- seq(from=min(grdVoorst$s1)+w, to=min(grdVoorst$s1)+5*w, by=w)+cellsize/2
grdVoorst$zone <- findInterval(grdVoorst$s1, s1bnd)

#construct all clusters
s1local <- grdVoorst$s1-min(grdVoorst$s1)
s2local <- grdVoorst$s2-min(grdVoorst$s2)
#set spacing of points within clusters
spacing <- 100 
s1f <- as.factor(s1local%%spacing)
s2f <- as.factor(s2local)
grdVoorst$cluster <- as.character(interaction(s1f,s2f, as.factor(grdVoorst$zone)))

tz_pop <- sum(grdVoorst$z)
tz_cl <- tapply(grdVoorst$z, INDEX=grdVoorst$cluster, FUN=sum)

M_cl <- tapply(grdVoorst$z, INDEX=grdVoorst$cluster, FUN=length)
p <- M_cl/sum(M_cl)
n <- 6

mz_pop <-mean(grdVoorst$z)
mz_cl <- tapply(grdVoorst$z, INDEX=grdVoorst$cluster, FUN=mean)
print(v_mz_cl_true <- sum(p*(mz_cl-mz_pop)^2)/n)

m_n <- n*sum(p*M_cl)
print(v_mz_SI_true <- var(grdVoorst$z)/m_n)
