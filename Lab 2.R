rm(list=ls())

library(PrevMap)
library(lme4)

chikwawa <- read.csv("lab_data.csv")

attach(chikwawa) # attach variables to use

# Extraction of the elevation data
elev.data <- unique(chikwawa[,c("web_x","web_y","elevation")]) 

# Testing correlation 
check.corr <- spat.corr.diagnostic(log(elevation) ~ 1, coords = ~web_x+web_y,
                                   data=elev.data,
                                   likelihood = "Gaussian",
                                   lse.variogram = TRUE,
                                   uvec = seq(2,20,length=15),
                                   which.test = "variogram")

# Parameter estimation
fit.mle <- 
linear.model.MLE(log(elevation) ~ 1, coords=~web_x+web_y,
                 data=elev.data,
                 start.cov.pars = 5.88,
                 fixed.rel.nugget = 0,
                 kappa=0.5) # this indicated the fitting of an exponential correlation
                 
summary(fit.mle,log.cov.pars=FALSE)


# Validation of the correlation function
lm.diag <- variog.diagnostic.lm(fit.mle,n.sim=10000,which.test = "variogram")

sigma2.hat <- coef(fit.mle)["sigma^2"]
phi.hat <- coef(fit.mle)["phi"]
curve(sigma2.hat*(1-exp(-x/phi.hat)),add=TRUE,col=2)

