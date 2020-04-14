# ---------------------------------------- #
# Cluster Detection of Cumulative SARS-CoV-2 incidence rate per 100,000 in Georgia, U.S.A.
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

library(DCluster)
library(sp)
library(spdep)

####################
# Data Importation #
####################

load(file = "data/CoV_GA.Rdata")

######################
# Spatial Clustering #
######################

# Local Moran's I
Y <- CoV_GA_UTM17N@data$cumrate
nb <- spdep::poly2nb(CoV_GA_UTM17N) # queen = TRUE
col.W <- spdep::nb2listw(nb, style = "B")

I.local <- spdep::localmoran(Y, col.W)
CoV_GA_UTM17N$I.local <- I.local[,1]
CoV_GA_UTM17N$I.local.p <- I.local[,5]
CoV_GA_UTM17N$I.local.p_bonf <- p.adjust(I.local[,5], method = "bonferroni")
CoV_GA_UTM17N$I.local.p_holm <- p.adjust(I.local[,5], method = "holm")
CoV_GA_UTM17N$I.local.p_fdr <- p.adjust(I.local[,5], method = "fdr")

png(file = "figures/Georgia_LocalICumulativeRate_stat.png", height = 1000*f, width = 1000*f)
sp::spplot(CoV_GA_UTM17N, 
           "I.local",
           main = "Local Moran's I statistic",
           col.regions = grey.colors(256),
           par.settings = list(axis.line = list(col =  'transparent'))
           )
dev.off()

CoV_GA_UTM17N$plot_junk <- cut(CoV_GA_UTM17N$I.local.p, c(0,0.01,0.05,0.1,1))
png(file = "figures/Georgia_LocalICumulativeRate_pvalues.png", height = 1000*f, width = 1000*f)
sp::spplot(CoV_GA_UTM17N, 
           "plot_junk", 
           main = "Local Moran's I p-value: Uncorrected",
           col.regions = rev(grey.colors(256)),
           par.settings = list(axis.line = list(col =  'transparent'))
           )
dev.off()

CoV_GA_UTM17N$plot_junk <- cut(CoV_GA_UTM17N$I.local.p_bonf, c(0,0.01,0.05,0.1,1))
png(file = "figures/Georgia_LocalICumulativeRate_bonferroni.png", height = 1000*f, width = 1000*f)
sp::spplot(CoV_GA_UTM17N, 
           "plot_junk", 
           main = "Local Moran's I p-value: Bonferroni Corrected",
           col.regions = rev(grey.colors(256)),
           par.settings = list(axis.line = list(col =  'transparent'))
           )
dev.off()

CoV_GA_UTM17N$plot_junk <- cut(CoV_GA_UTM17N$I.local.p_holm, c(0,0.01,0.05,0.1,1))
png(file = "figures/Georgia_LocalICumulativeRate_holm.png", height = 1000*f, width = 1000*f)
sp::spplot(CoV_GA_UTM17N, 
           "plot_junk", 
           main = "Local Moran's I p-value: Holm Adjustment",
           col.regions = rev(grey.colors(256)),
           par.settings = list(axis.line = list(col =  'transparent'))
           )
dev.off()

CoV_GA_UTM17N$plot_junk <- cut(CoV_GA_UTM17N$I.local.p_fdr, c(0,0.01,0.05,0.1,1))
png(file = "figures/Georgia_LocalICumulativeRate_fdr.png", height = 1000*f, width = 1000*f)
sp::spplot(CoV_GA_UTM17N, 
           "plot_junk", 
           main = "Local Moran's I p-value: False Discovery Rate",
           col.regions = rev(grey.colors(256)),
           par.settings = list(axis.line = list(col =  'transparent'))
           )
dev.off()

# Permutation test of local Moran's I
n.iter <- 1000
I.keep <- matrix(NA, ncol = n.iter, nrow = length(Y))
for(i in 1:n.iter){
  if (i %% 200 == 0){print(i)}
  I.keep[,i] <- spdep::localmoran(sample(Y, length(Y), replace = T), 
                                  col.W)[,1]
  }

par(mfrow = c(1,2))
qqnorm(I.keep[1,]); qqline(I.keep[1,])
qqnorm(I.keep[20,]); qqline(I.keep[20,])
par(op)

I.obs <- spdep::localmoran(Y, col.W)[,1]
P.perm <- apply(sweep(I.keep, 1, I.obs, ">"), 1, mean)
P.perm[P.perm == 0] <- 1/n.iter
CoV_GA_UTM17N$plot_junk <- cut(p.adjust(P.perm, "holm"), c(0.01, 0.05, 0.1, 1))

png(file = "figures/Georgia_permLocalICumulativeRate.png", height = 1000*f, width = 1000*f)
sp::spplot(CoV_GA_UTM17N, 
           "plot_junk", 
           main = "Local Moran's I p-value: Permutation with Holm Adjustment",
           col.regions = rev(grey.colors(256)),
           par.settings = list(axis.line = list(col =  'transparent'))
           )
dev.off()

#####################################
# Disease Mapping Cluster Detection #
#####################################

dismap <- data.frame(Observed = CoV_GA_UTM17N$cumulative,
                     Pop = CoV_GA_UTM17N$Pop2010
                     )
theta_0 <- sum(dismap$Observed)/sum(dismap$Pop)
dismap$Expected <- dismap$Pop * theta_0

# Chi-Squared Test
## Asymptotic test
DCluster::achisq.stat(dismap) # significant 
## Simulation-based
### Poisson
DCluster::achisq.test(Observed ~ offset(log(Expected)),
                      data = dismap,
                      model = "poisson",
                      R = 1000
                      ) # significant 
### Multinomial
DCluster::achisq.test(Observed ~ offset(log(Expected)),
                      data = dismap,
                      model = "multinom",
                      R = 1000
                      ) # significant 
### Negative Binomial
DCluster::achisq.test(Observed ~ offset(log(Expected)),
                      data = dismap,
                      model = "negbin",
                      R = 1000
                      ) # not significant 

# Potthoff-Whittinghill Test
## Asymptotic test
DCluster::pottwhitt.stat(dismap) # significant 
## Simulation-based
### Poisson
DCluster::pottwhitt.test(Observed ~ offset(log(Expected)),
                         data = dismap,
                         model = "poisson",
                         R = 1000
                         ) # significant 
### Multinomial
DCluster::pottwhitt.test(Observed ~ offset(log(Expected)),
                         data = dismap,
                         model = "multinom",
                         R = 1000
                         ) # significant 
### Negative Binomial
DCluster::pottwhitt.test(Observed ~ offset(log(Expected)),
                         data = dismap,
                         model = "negbin",
                         R = 1000
                         ) # not significant 

# Moran's I
nb <- spdep::poly2nb(CoV_GA_UTM17N)
col.W <- spdep::nb2listw(nb)
## Asymptotic test
DCluster::moranI.stat(dismap, 
                      listw = col.W,
                      n = nrow(dismap),
                      S0 = spdep::Szero(col.W)
                      )
## Simulation-based
### Poisson
DCluster::moranI.test(Observed ~ offset(log(Expected)),
                      data = dismap,
                      model = "poisson",
                      R = 1000,
                      listw = col.W,
                      n = nrow(dismap),
                      S0 = spdep::Szero(col.W)
                      ) # significant 
### Multinomial
### Poisson
DCluster::moranI.test(Observed ~ offset(log(Expected)),
                      data = dismap,
                      model = "multinom",
                      R = 1000,
                      listw = col.W,
                      n = nrow(dismap),
                      S0 = spdep::Szero(col.W)
                      ) # significant 
### Negative Binomial
DCluster::moranI.test(Observed ~ offset(log(Expected)),
                      data = dismap,
                      model = "negbin",
                      R = 1000,
                      listw = col.W,
                      n = nrow(dismap),
                      S0 = spdep::Szero(col.W)
                      ) # significant 

# Tango's T
coords <- sp::coordinates(CoV_GA_UTM17N)
dlist <- spdep::dnearneigh(coords, 0, Inf)
dlist <- spdep::include.self(dlist)
dlist.d <- spdep::nbdists(dlist, coords, longlat = F)
for (i in 1:length(dlist.d)){dlist.d[[i]][i] <- 0}
phi <- 20
col.W.tango <- spdep::nb2listw(dlist, 
                               glist = lapply(dlist.d, function(x){exp(-x/phi)}),
                               style = "C"
                               )
## Asymptotic
DCluster::tango.stat(dismap, listw = col.W.tango)
## Simulation-based
### Poisson
DCluster::tango.test(Observed ~ offset(log(Expected)),
                     data = dismap,
                     model = "poisson",
                     R = 1000,
                     listw = col.W.tango
                     ) # significant 
### Multinomial
DCluster::tango.test(Observed ~ offset(log(Expected)),
                     data = dismap,
                     model = "multinom",
                     R = 1000,
                     listw = col.W.tango
                     ) # significant 
### Negative Binomial
DCluster::tango.test(Observed ~ offset(log(Expected)),
                     data = dismap,
                     model = "negbin",
                     R = 1000,
                     listw = col.W.tango
                     ) # insignificant 

# Local Clusters
## Convert to kilometers (must use projected UTM)
dismap$x <- sp::coordinates(CoV_GA_UTM17N)[,1]/1000
dismap$y <- sp::coordinates(CoV_GA_UTM17N)[,2]/1000

## estimate MLE parameters for the negative binomial model
mle <- DCluster::calculate.mle(dismap, model = "negbin")

## Openshaw's GAM
opgam.fit <- DCluster::opgam(data = dismap,
                             radius = 50,
                             alpha = 0.01,
                             iscluster = opgam.iscluster.negbin,
                             mle = mle
                             )

png(file = "figures/Georgia_OpenshawCumulativeRate1.png", height = 1000*f, width = 1000*f)
plot(CoV_GA_UTM17N, main = "Openshaw's GAM of Cumulative COVID-19 Rate per 100,000 (radius = 50 km; alpha = 0.01)")
lines(opgam.fit[,1]*1000,
      opgam.fit[,2]*1000,
      type = "p",
      pch = "x",
      cex = 0.8,
      col = 2
      )
dev.off()

opgam.fit <- DCluster::opgam(data = dismap,
                             radius = 75,
                             alpha = 0.05,
                             iscluster = opgam.iscluster.negbin,
                             mle = mle
                             )

png(file = "figures/Georgia_OpenshawCumulativeRate2.png", height = 1000*f, width = 1000*f)
plot(CoV_GA_UTM17N, main = "Openshaw's GAM of Cumulative COVID-19 Rate per 100,000 (radius = 75 km; alpha = 0.05)")
lines(opgam.fit[,1]*1000,
      opgam.fit[,2]*1000,
      type = "p",
      pch = "x",
      cex = 0.8,
      col = 2
      )
dev.off()

## Kulldorf and Nagarwalle
mle <- DCluster::calculate.mle(dismap, model = "negbin")
grid.search <- sp::coordinates(CoV_GA_UTM17N)/1000

opgam.fit <- DCluster::opgam(data = dismap,
                             thegrid = grid.search,
                             alpha = 0.05,
                             iscluster = kn.iscluster,
                             model = "negbin",
                             fractpop = 0.1,
                             R = 500,
                             mle = mle
                             )

png(file = "figures/Georgia_KNCumulativeRate1.png", height = 1000*f, width = 1000*f)
plot(CoV_GA_UTM17N, main = "Kulldorf and Nagarwalle of Cumulative COVID-19 Rate per 100,000 (fraction of population = 0.1; alpha = 0.05)")
lines(opgam.fit[,1]*1000,
      opgam.fit[,2]*1000,
      type = "p",
      pch = 16,
      cex = 1.2,
      col = 2
      )
dev.off()

opgam.fit <- DCluster::opgam(data = dismap,
                             thegrid = grid.search,
                             alpha = 0.05,
                             iscluster = kn.iscluster,
                             model = "negbin",
                             fractpop = 0.05,
                             R = 500,
                             mle = mle
                             )

png(file = "figures/Georgia_KNCumulativeRate2.png", height = 1000*f, width = 1000*f)
plot(CoV_GA_UTM17N, main = "Kulldorf and Nagarwalle of Cumulative COVID-19 Rate per 100,000 (fraction of population = 0.05; alpha = 0.05)")
lines(opgam.fit[,1]*1000,
      opgam.fit[,2]*1000,
      type = "p",
      pch = 16,
      cex = 1.2,
      col = 2
)
dev.off()

## Besag and Newell

opgam.fit <- DCluster::opgam(data = dismap,
                             thegrid = grid.search,
                             alpha = 0.05,
                             iscluster = bn.iscluster,
                             model = "negbin",
                             k = 500,
                             R = 500,
                             mle = mle
                             )

png(file = "figures/Georgia_BNCumulativeRate1.png", height = 1000*f, width = 1000*f)
plot(CoV_GA_UTM17N, main = "Besag and Newell of Cumulative COVID-19 Rate per 100,000 (size of cluster = 500; alpha = 0.05)")
lines(opgam.fit[,1]*1000,
      opgam.fit[,2]*1000,
      type = "p",
      pch = 16,
      cex = 1.2,
      col = 2
      )
dev.off()

opgam.fit <- DCluster::opgam(data = dismap,
                             thegrid = grid.search,
                             alpha = 0.05,
                             iscluster = bn.iscluster,
                             model = "negbin",
                             k = 2000,
                             R = 500,
                             mle = mle
                             )

png(file = "figures/Georgia_BNCumulativeRate2.png", height = 1000*f, width = 1000*f)
plot(CoV_GA_UTM17N, main = "Besag and Newell of Cumulative COVID-19 Rate per 100,000 (size of cluster = 2,000; alpha = 0.05)")
lines(opgam.fit[,1]*1000,
      opgam.fit[,2]*1000,
      type = "p",
      pch = 16,
      cex = 1.2,
      col = 2
      )
dev.off()
# -------------------- END OF CODE -------------------- #