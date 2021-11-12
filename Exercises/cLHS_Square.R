library(clhs)
library(ggplot2)

#Load simulated square
load(file="../data/SimulatedSquare.RData")

#Use covariate $x$ together with the spatial coordinates s1 and s2 to select a conditioned Latin hypercube sample of 16 points

set.seed(314)
res <- clhs(grid[, c(1,2,3)], size = 16, tdecrease = 0.95, iter = 10000, progress = FALSE, simple = FALSE)
index <- res$index_samples
myclhsample <- grid[index, ]

ggplot(data=grid) +
  geom_tile(mapping = aes(x = s1, y = s2, fill = x))+ 
  geom_tile(data=myclhsample,mapping = aes(x = s1, y = s2), width=1,height=1,size=1,colour="red",fill = NA) +
  scale_x_continuous(name = "Easting") +
  scale_y_continuous(name = "Northing") + 
  scale_fill_continuous(name="x",type="viridis") +
  coord_fixed()

probs <- seq(from=0,to=1,length.out = 16 + 1)
breaks <- apply(grid[,c(1,2,3)],MARGIN=2,FUN=function(x) quantile(x,probs=probs))
counts <- lapply(1:3, function (i) 
  hist(myclhsample[, i], breaks[,i], plot = FALSE)$counts
)

countslf <- data.frame(counts=unlist(counts))
countslf$covariate <- rep(c("s1","s2","x"),each=16)
countslf$stratum <- rep(1:16,times=3)
ggplot(countslf) +
  geom_point(mapping = aes(x=stratum,y = counts), colour = "black",size=1) +
  facet_wrap(~covariate) +
  scale_x_continuous(name = "Stratum") +
  scale_y_continuous(name = "Sample size",breaks=c(0,1,2))

trace <- res$obj
tracedf <- data.frame(trace=trace)
ggplot(tracedf) +
  geom_line(mapping = aes(x=1:nrow(tracedf),y = trace), colour = "black",size=0.8) +
  scale_x_continuous(name = "Iteration") +
  scale_y_continuous(name = "Criterion")
tail(trace, n=1)
