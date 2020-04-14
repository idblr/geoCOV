# ---------------------------------------- #
# Spatial Lag Regression of Cumulative SARS-CoV-2 incidence rate per 100,000 in Georgia, U.S.A.
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

####################
# Data Importation #
####################

load(file = "data/CoV_GA.Rdata")

##########################
# Spatial Lag Regression #
##########################

Y <- CoV_GA_UTM17N@data$cumrate # cumulative incidence rate per 100,000
W <- spdep::nb2mat(nb, style = "W")
WY <- W %*% Y

png(file = "figures/Georgia_SpatiallyLaggedCumulativeRate.png", height = 1000*f, width = 1000*f)
par(pty = "s")
plot(WY ~ Y, 
     ylab = "Spatially Lagged Y",
     xlab = "Cumulative Rate per 100,000 (Y)",
     cex.lab = 1,
     cex.axis = 1,
     cex = 1,
     col = 4
)
abline(h = mean(WY), lty = 2)
abline(v = mean(Y), lty = 2)
abline(0,1)
dev.off()

# another method
spdep::moran.plot(Y, spdep::nb2listw(nb))

# Calculate Moran's I by hand and perform asymptotic hypothesis test
ybar <- mean(Y)
r <- Y - ybar
I <- sum(r %*% t(r) %*% W) / sum(r^2) * nrow(W) / sum(W)
col.W <- spdep::nb2listw(nb, style = "W")

# default performs non-Gaussian test under randomization
spdep::moran.test(Y, col.W) # significant 0.57

# assume Y is Gaussian
spdep::moran.test(Y, col.W, randomisation = F) # significant 0.57

# Moran's I hypothesis test by permutation
I.perm <- spdep::moran.mc(Y, col.W, 10000)
I.perm # significant 0.57

png(file = "figures/Georgia_MoransICumulativeRate.png", height = 1000*f, width = 1000*f)
par(pty = "s")
hist(I.perm$res,
     main = "Moran's I under Null Hypothesis",
     xlab = "",
     ylab = "",
     xlim = c(-0.2,0.6),
     ylim = c(0,5000)
)
abline(v = 0.570, col = 2, lwd = 4)
text(0.5, 4000, "Observed \n Moran's I", col = 2)
dev.off()

# Geary's C
spdep::geary.test(Y, col.W) # significant 0.456

# Geary's C hypothesis test by permutation
C.perm <- spdep::geary.mc(Y, col.W, 10000)
png(file = "figures/Georgia_GearyCCumulativeRate.png", height = 1000*f, width = 1000*f)
par(pty = "s")
hist(C.perm$res,
     main = "Geary's C under Null Hypothesis",
     xlab = "",
     ylab = "",
     xlim = c(0.4, 1.4),
     ylim = c(0, 3000)
)
abline(v = 0.456, col = 2, lwd = 4)
text(0.5, 2500, "Observed \n Geary's C", col = 2)
dev.off()
# -------------------- END OF CODE -------------------- #