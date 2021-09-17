load(file="data/SpainPortugal_TAS.RData")


n <- 100 #size of simple random samples per time
m <- 50 #size of matched sample

#compute population variance-covariance matrix of data
S2_z <- var(df[,-c(1,2)])

#compute variance of change of mean using pi estimators to estimate the spatial means of 2004 and 2019
(se_d_HT <- sqrt(S2_z[1,1]/n + S2_z[4,4]/n - 2*(m*S2_z[1,4]/n^2)))



#Supplemented panel
t <- 1:4
C_SS <- S2_z/(n/2)
C_IS <- C_SS
C_IS[row(C_IS)!=col(C_IS)] <- 0
zeroes <- matrix(0, nrow=4, ncol=4)
C <- rbind(cbind(C_SS,zeroes),cbind(zeroes,C_IS))
X <- rbind(diag(4),diag(4))
Cinv <- solve(C)
XCXinv <- solve(crossprod(X, Cinv)%*%X)
w <- c(-1,0,0,1)
(se_d_GLS <- sqrt(t(w)%*%XCXinv%*%w))
