# ---------------------------------------- #
# Disease mapping of SARS-CoV-2 incidence rate in Georgia, U.S.A.
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

library(colorRamps)
library(sp)
library(spdep)
library(CARBayes)

###################
# Custom Function #
###################
line01 <- function(x, y, ...){
  abline(0, 1, col = 4, lwd = 2); points(x, y, col = 2)
  }

####################
# Data Importation #
####################

load(file = "data/CoV_GA.Rdata")

###################
# Disease Mapping #
###################

Y <- CoV_GA_UTM17N$cumulative
P <- CoV_GA_UTM17N$Pop2010
X <- CoV_GA_UTM17N$pernoncauc

nb <- spdep::poly2nb(CoV_GA_UTM17N)
W <- spdep::nb2mat(nb, style = "B") # binary

# MLE Estimate
u.mle <- log(Y / P) # incidence rate
se.mle <- sqrt(1 / Y) # standard error

# Exchangeable Model
fit.exch <- CARBayes::S.CARbym(Y ~ offset(log(P)),
                               family = "poisson",
                               W = W,
                               prior.tau2 = c(1e50, 1), # uniformed
                               burnin = 10000,
                               n.sample = 100000,
                               verbose = F
                               )
names(fit.exch) # see values
## Trace plots
plot(fit.exch$samples[["beta"]])
plot(fit.exch$samples[["sigma2"]])
## Summary
fit.exch$summary.results
## Extract the posterior sample
### Also add in the intercept so they are comparable to MLE estimates
post.exch <- sweep(fit.exch$samples[["psi"]], 1, fit.exch$samples[["beta"]], "+")

# CAR-only model
fit.car <- CARBayes::S.CARbym(Y ~ offset(log(P)),
                              family = "poisson",
                              W = W,
                              prior.sigma2 =  c(1e50, 1),
                              burnin = 10000,
                              n.sample = 100000,
                              verbose = F
                              )
post.car <- sweep(fit.car$samples[["psi"]], 1, fit.car$samples[["beta"]], "+")

# Convolution BYM. Use the default priors
fit.conv <- CARBayes::S.CARbym(Y ~ offset(log(P)),
                               family = "poisson",
                               W = W,
                               burnin = 10000,
                               n.sample = 100000,
                               verbose = F
                               )
post.conv <- sweep(fit.conv$samples[["psi"]], 1, fit.conv$samples[["beta"]], "+")

# Save results
CoV_GA_UTM17N@data$rr.mle <- u.mle
CoV_GA_UTM17N@data$se.mle <- se.mle
CoV_GA_UTM17N@data$rr.exch <- colMeans(post.exch)
CoV_GA_UTM17N@data$se.exch <- apply(post.exch, 2, sd)
CoV_GA_UTM17N@data$rr.car <- colMeans(post.car)
CoV_GA_UTM17N@data$se.car <- apply(post.car, 2, sd)
CoV_GA_UTM17N@data$rr.conv <- colMeans(post.conv)
CoV_GA_UTM17N@data$se.conv <- apply(post.conv, 2, sd)

pairs(CoV_GA_UTM17N@data[c("rr.mle", "rr.exch", "rr.car", "rr.conv")], panel = line01)
pairs(CoV_GA_UTM17N@data[c("se.mle", "se.exch", "se.car", "se.conv")], panel = line01)

changeRR <- CoV_GA_UTM17N@data$rr.mle - CoV_GA_UTM17N@data$rr.conv
changeSE <- CoV_GA_UTM17N@data$se.mle - CoV_GA_UTM17N@data$se.conv

par (mfrow = c(1,2), mar = c(4, 4, 1, 1))
plot(changeRR ~ sqrt(1/Y), 
     ylab = "Change in log RR",
     xlab = "SQRT (1/Y)"); abline(h = 0, col = 4)
plot(changeRR ~ log(P), 
     ylab = "Change in log RR",
     xlab = "log Population"); abline(h = 0, col = 4)

par (mfrow = c(1,2), mar = c(4, 4, 1, 1))
plot(changeSE ~ sqrt(1/Y), 
     ylab = "Change in SE of log RR",
     xlab = "SQRT (1/Y)"); abline(h = 0, col = 4)
plot(changeSE ~ log(P), 
     ylab = "Change in SE of log RR",
     xlab = "log Population"); abline(h = 0, col = 4)

par (mfrow = c(1,2), mar = c(4, 4, 1, 1))
plot(as.numeric(fit.conv$samples[["tau2"]]),
     type = "l",
     ylab = "Tau^2"
     )
plot(as.numeric(fit.conv$samples[["sigma2"]]),
     type = "l",
     ylab = "Sigma^2"
     )

png(file = "figures/Georgia_CAR_rr.png", height = 1000*f, width = 1000*f)
sp::spplot(CoV_GA_UTM17N, 
           z = c("rr.mle", "rr.exch", "rr.car", "rr.conv"),
           main = "Log Incidence",
           col.regions = colorRamps::matlab.like(256),
           par.settings = list(axis.line = list(col =  'transparent'))
)
dev.off()

png(file = "figures/Georgia_CAR_se.png", height = 1000*f, width = 1000*f)
sp::spplot(CoV_GA_UTM17N, 
           z = c("se.mle", "se.exch", "se.car", "se.conv"),
           main = "Standard Deviation Log Incidence",
           col.regions = colorRamps::matlab.like(256),
           par.settings = list(axis.line = list(col =  'transparent'))
)
dev.off()

# BYM with percent non-caucasian
fit.noncauc <- CARBayes::S.CARbym(Y ~ offset(log(P)) + X,
                                  family = "poisson",
                                  W = W,
                                  burnin = 10000,
                                  n.sample = 200000,
                                  thin = 5, 
                                  verbose = F
)
post.noncauc <- fit.noncauc$samples[["psi"]]
CoV_GA_UTM17N@data$rr.noncauc <- colMeans(post.noncauc)
CoV_GA_UTM17N@data$se.noncauc <- apply(post.noncauc, 2, sd)
CoV_GA_UTM17N@data$p.higher.noncauc <- apply(post.noncauc > 0, 2, mean)
CoV_GA_UTM17N@data$p.lower.noncauc <- apply(post.noncauc < 0, 2, mean)

# fit standard poisson and quasi-poisson models
fit.cauc.pois <- glm(Y ~ offset(log(P)) + X, family = "poisson"); summary(fit.cauc.pois) # significant
fit.cauc.qpois <- glm(Y ~ offset(log(P)) + X, family = "quasipoisson"); summary(fit.cauc.qpois) # significant

png(file = "figures/Georgia_CAR_noncauc_rrr.png", height = 1000*f, width = 1000*f)
sp::spplot(CoV_GA_UTM17N, 
           z = c("rr.noncauc", "se.noncauc"),
           main = "Residual Random Effects",
           col.regions = colorRamps::matlab.like(256),
           par.settings = list(axis.line = list(col =  'transparent'))
)
dev.off()

png(file = "figures/Georgia_CAR_noncauc_prob0.png", height = 1000*f, width = 1000*f)
sp::spplot(CoV_GA_UTM17N, 
           z = c("p.higher.noncauc", "p.lower.noncauc"),
           main = "Probability of > 0",
           col.regions = colorRamps::matlab.like(256),
           par.settings = list(axis.line = list(col =  'transparent'))
)
dev.off()
# -------------------- END OF CODE -------------------- #