library(clhs)
library(ggplot2)

#load grid with covariates
grdHunterValley <- readRDS(file = "../data/grdHunterValley.rds")
names(grdHunterValley)

n <- 50
set.seed(314)
res <- clhs(grdHunterValley[,c(3,4,5,6,7)], size=n, tdecrease=0.8, length.cycle=100, progress=TRUE, simple=FALSE)

str(res)
plot(res$obj)

mysample <- res$sampled_data
mysample$s1 <- grdHunterValley$s1[res$index_samples]
mysample$s2 <- grdHunterValley$s2[res$index_samples]

#Plot the selected points on top of one of the covariates
ggplot(data=grdHunterValley) +
  geom_tile(mapping = aes(x = s1, y = s2, fill = cti))+  
  geom_point(data = mysample, mapping = aes(x = s1, y = s2), colour = "red") +
  scale_x_continuous(name = "Easting (km)") +
  scale_y_continuous(name = "Northing (km)") +    
  scale_fill_viridis_c(name="cti")+
  coord_fixed()

#Make scatter plot
ggplot(data=grdHunterValley) +
  geom_point(mapping = aes(x = cti, y = elevation_m), colour = "black",size=1, alpha=0.5) +
  geom_point(data=mysample, mapping = aes(x = cti, y = elevation_m), colour = "red", size=2) +
  scale_x_continuous(name = "cti") +
  scale_y_continuous(name = "elevation")