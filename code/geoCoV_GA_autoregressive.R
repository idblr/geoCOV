# ---------------------------------------- #
# Autoregressive models of Cumulative SARS-CoV-2 incidence rate per 100,000 in Georgia, U.S.A.
# ---------------------------------------- #
#
# Created by: Ian Buller, Ph.D., M.A. (GitHub: @idblr)
# Created on: April 13, 2020
#
# Recently modified by:
# Recently modified on:
#
# Notes:
# A) 04/13/2020 (IB) - Builds upon data generated in the "geoCoV_CA.R" file
# B) 04/13/2020 (IB) - Credit to Dr. Howard Chang (BIOS737 at Emory University)
# ---------------------------------------- #

############
# Packages #
############

library(sp)
library(spdep)
library(spatialreg)

####################
# Data Importation #
####################

load(file = "data/CoV_GA.Rdata")

#########################
# Autoregressive Models #
#########################

## Simultaneous Autoregressive Models
nb <- spdep::poly2nb(CoV_GA_UTM17N)
nb <- spdep::nblag(nb,2)
W <- spdep::nb2listw(nb[[1]], style = "W")
W2 <- spdep::nb2listw(nb[[2]], style = "W")

Y <- CoV_GA_UTM17N$cumrate
X <- CoV_GA_UTM17N$pernoncauc
P <- CoV_GA_UTM17N$Pop2010

### Standard regression with and wihtout population weights
fit0 <- lm(Y ~ X); summary(fit0)
fit0.p <- lm(Y ~ X, weight = P); summary(fit0.p)

### Spatial Error Models
fit1 <- spatialreg::spautolm(Y ~ X, listw = W); summary(fit1)
fit2 <- spatialreg::spautolm(Y ~ X, listw = W2); summary(fit2)
fit1.p <- spatialreg::spautolm(Y ~ X, listw = W, weights = P); summary(fit1.p)
fit2.p <- spatialreg::spautolm(Y ~ X, listw = W2, weights = P); summary(fit2.p)

### Spatial Lag Model
fit1.lag <- spatialreg::lagsarlm(Y ~ X, listw = W); summary(fit1.lag)
fit2.lag <- spatialreg::lagsarlm(Y ~ X, listw = W2); summary(fit2.lag)
fit1.durbin <- spatialreg::lagsarlm(Y ~ X, listw = W, type = "mixed"); summary(fit1.durbin)
fit2.durbin <- spatialreg::lagsarlm(Y ~ X, listw = W2, type = "mixed"); summary(fit2.durbin)

## Conditional Autoregressive Models
nb <- spdep::poly2nb(CoV_GA_UTM17N)
nb <- spdep::nblag(nb,2)
W <- spdep::nb2listw(nb[[1]], style = "B")
W2 <- spdep::nb2listw(nb[[2]], style = "B")

fit1.car <- spatialreg::spautolm(Y ~ X, family = "CAR", listw = W); summary(fit1.car)
fit2.car <- spatialreg::spautolm(Y ~ X, family = "CAR", listw = W2); summary(fit2.car)

# Examin impacts of parameter choice on SAR and CAR
nb <- spdep::poly2nb(CoV_GA_UTM17N)
W <- spdep::nb2mat(nb, style = "B") # binary
W.row <- spdep::nb2mat(nb, style = "W") # row-standardized
n <- nrow(W)
D <- diag(apply(W, 1, sum))

lambdas <- seq(0, 0.99, by = 0.01)
cov.SAR <- cov.CAR <- array(NA, c(n,n,length(lambdas)))
cor.1st.SAR <- cor.1st.CAR <- matrix(NA, n, length(lambdas))
var.SAR <- var.CAR <- matrix(NA, n, length(lambdas))

for (i in 1:length(lambdas)) {
  lambda.i <- lambdas[i]
  cov.SAR.i <- solve(diag(n) - lambda.i * W.row) %*% solve(diag(n) - lambda.i * t(W.row))
  cov.CAR.i <- solve(D - lambda.i * W)
  cor.SAR.i <- cov2cor(cov.SAR.i) 
  cor.CAR.i <- cov2cor(cov.CAR.i) 
  cov.SAR[,,i] <- cov.SAR.i
  cov.CAR[,,i] <- cov.CAR.i
  cor.1st.SAR[,i] <- rowSums(cor.SAR.i * W) / diag(D)
  cor.1st.CAR[,i] <- rowSums(cor.CAR.i * W) / diag(D)
  var.SAR[,i] <- diag(cov.SAR.i)
  var.CAR[,i] <- diag(cov.CAR.i)
}

plot(0, 0, 
     ylim = c(0,1), xlim = c(0,1),
     xlab = "Spatial Parameter", ylab = "County-specific Correlation with 1st Neighbor",
     main = "SAR",
     type = "n")
for (i in 1:n) {lines(cor.1st.SAR[i,] ~ lambdas, col = rgb(0, 0, 1, 0.2))}

plot(0, 0, 
     ylim = c(0, 10), xlim = c(0, 1),
     xlab = "Spatial Parameter", ylab = "County-specific Standard Deviation",
     main = "SAR",
     type = "n")
for (i in 1:n) {lines(sqrt(var.SAR[i,]) ~ lambdas, col = rgb(0, 0, 1, 0.2))}

plot(0, 0, 
     ylim = c(0,1), xlim = c(0,1),
     xlab = "Spatial Parameter", ylab = "County-specific Correlation with 1st Neighbor",
     main = "CAR",
     type = "n")
for (i in 1:n) {lines(cor.1st.CAR[i,] ~ lambdas, col = rgb(0, 0, 1, 0.2))}

plot(0, 0, 
     ylim = c(0, 2), xlim = c(0, 1),
     xlab = "Spatial Parameter", ylab = "County-specific Standard Deviation",
     main = "CAR",
     type = "n")
for (i in 1:n) {lines(sqrt(var.CAR[i,]) ~ lambdas, col = rgb(0, 0, 1, 0.2))}
# -------------------- END OF CODE -------------------- #