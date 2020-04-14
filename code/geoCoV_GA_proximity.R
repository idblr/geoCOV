# ---------------------------------------- #
# Defining Spatial Proximity of Georgia, U.S.A.
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

##############################
# Defining Spatial Proximity #
##############################

nb <- spdep::poly2nb(CoV_GA_UTM17N) # queen = TRUE
nb2 <- spdep::poly2nb(CoV_GA_UTM17N, queen = FALSE)
summary(nb)

W <- spdep::nb2mat(nb, style = "B")
centroids <- sp::coordinates(CoV_GA_UTM17N)

f <- 1 # exansion factor
png(file = "figures/Georgia_Neighbors.png", height = 1000*f, width = 1000*f)
plot(CoV_GA_UTM17N)
plot(nb, centroids, add = T, col = "blue", lwd = 1*f)
plot(nb2, centroids, add = T, col = "red", lwd = 1*f)
title(main = "Blue lines = Queen's case, Red lines = Rook's case", cex = 1*f)
dev.off()

## Using kth-nearest distance
nb_k1 <- spdep::knn2nb(spdep::knearneigh(centroids, k = 1, longlat = FALSE),
                       row.names = row.names(centroids)
                       )
nb_k2 <- spdep::knn2nb(spdep::knearneigh(centroids, k = 2, longlat = FALSE),
                       row.names = row.names(centroids)
                       )

f <- 1 # exansion factor
png(file = "figures/Georgia_KthNearestNeighbors.png", height = 1000*f, width = 1000*f)
plot(CoV_GA_UTM17N)
plot(nb_k2, centroids, add = T, col = "blue", lwd = 1*f)
plot(nb_k1, centroids, add = T, col = "red", lwd = 1*f)
title(main = "Blue lines = two nearest neighbors, Red lines = one nearset neighbor", cex = 1*f)
dev.off()

## Using buffer distance
nb_b1 <- spdep::dnearneigh(centroids, d1 = 0, d2 = 30000, longlat = FALSE)
nb_b2 <- spdep::dnearneigh(centroids, d1 = 0, d2 = 60000, longlat = FALSE)

f <- 1 # exansion factor
png(file = "figures/Georgia_BufferNeighbors.png", height = 1000*f, width = 1000*f)
plot(CoV_GA_UTM17N)
plot(nb_b2, centroids, add = T, col = "blue", lwd = 1*f)
plot(nb_b1, centroids, add = T, col = "red", lwd = 1*f)
title(main = "Blue lines = 600 km apart, Red lines = 300 km apart", cex = 1*f)
dev.off()
# -------------------- END OF CODE -------------------- #