library(sswr)

n <- 100 #size of simple random samples per time
t <- 1:4 #sampling times (unit of time is five years)

#compute population variance-covariance matrix of data
S2_z <- var(grdSpainPortugal[, -c(1, 2)])

#Static-synchronous

C <- S2_z / n
#current mean
se_mz_cur_SS <- sqrt(C[4, 4])
#change of mean
w <- c(-1, 0, 0, 1)
se_d_mz_SS <- sqrt(t(w) %*% C %*% w)
#trend of mean
w <- (t - mean(t)) / sum((t - mean(t))^2)
se_mz_trend_SS <- sqrt(t(w) %*% C %*% w)
#space-time mean
w <- rep(1 / 4, 4)
se_mz_st_SS <- sqrt(t(w) %*% C %*% w)

#Independent-synchronous

C <- S2_z / n
C[row(C) != col(C)] <- 0
#current mean
se_mz_cur_IS <- sqrt(C[4, 4])
#change of mean
w <- c(-1, 0, 0, 1)
se_d_mz_IS <- sqrt(t(w) %*% C %*% w)
#trend of mean
w <- (t - mean(t)) / sum((t - mean(t))^2)
se_mz_trend_IS <- sqrt(t(w) %*% C %*% w)
#space-time mean
w <- rep(1 / 4, 4)
se_mz_st_IS <- sqrt(t(w) %*% C %*% w)

#Serially alternating

C <- S2_z / n
odd <- c(1, 3)
C[row(C) %in% odd & !(col(C) %in% odd)] <- 0
C[!(row(C) %in% odd) & col(C) %in% odd] <- 0

#current mean
se_mz_cur_SA <- sqrt(C[4, 4])
#change of mean
w <- c(-1, 0, 0, 1)
se_d_mz_SA <- sqrt(t(w) %*% C %*% w)
#trend of mean
w <- (t - mean(t)) / sum((t - mean(t))^2)
se_mz_trend_SA <- sqrt(t(w) %*% C %*% w)
#space-time mean
w <- rep(1 / 4, 4)
se_mz_st_SA <- sqrt(t(w) %*% C %*% w)

#Supplemented panel

C_SS <- S2_z / (n / 2)
C_IS <- C_SS
C_IS[row(C_IS) != col(C_IS)] <- 0
zeroes <- matrix(0, nrow = 4, ncol = 4)
C <- rbind(cbind(C_SS, zeroes), cbind(zeroes, C_IS))
X <- rbind(diag(4), diag(4))
XCXinv <- solve(crossprod(X, solve(C, X)))

#current mean
se_mz_cur_SP <- sqrt(XCXinv[4, 4])
#change of mean
w <- c(-1, 0, 0, 1)
se_d_mz_SP <- sqrt(t(w) %*% XCXinv %*% w)
#trend of mean
w <- (t - mean(t)) / sum((t - mean(t))^2)
se_mz_trend_SP <- sqrt(t(w) %*% XCXinv %*% w)
#space-time mean
w <- rep(1 / 4, 4)
se_mz_st_SP <- sqrt(t(w) %*% XCXinv %*% w)

#Rotating panel
C_b <- S2_z[1:2, 1:2] / (n / 2)
C_c <- S2_z[2:3, 2:3] / (n / 2)
C_d <- S2_z[3:4, 3:4] / (n / 2)
C <- matrix(data = 0, ncol = 8, nrow = 8)
C[2:3, 2:3] <- C_b
C[4:5, 4:5] <- C_c
C[6:7, 6:7] <- C_d
C[1, 1] <- S2_z[1, 1] / (n / 2)
C[8, 8] <- S2_z[4, 4] / (n / 2)
X <- matrix(c(rep(c(1, 0, 0, 0), 2),
              rep(c(0, 1, 0, 0), 2),
              rep(c(0, 0, 1, 0), 2),
              rep(c(0, 0, 0, 1), 2)), nrow = 8, ncol = 4, byrow = TRUE)
XCXinv <- solve(crossprod(X, solve(C, X)))

#current mean
se_mz_cur_RP <- sqrt(XCXinv[4, 4])
#change of mean
w <- c(-1, 0, 0, 1)
se_d_mz_RP <- sqrt(t(w) %*% XCXinv %*% w)
#trend of mean
w <- (t - mean(t)) / sum((t - mean(t))^2)
se_mz_trend_RP <- sqrt(t(w) %*% XCXinv %*% w)
#space-time mean
w <- rep(1 / 4, 4)
se_mz_st_RP <- sqrt(t(w) %*% XCXinv %*% w)
