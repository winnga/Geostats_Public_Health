rm(list=ls())

library(PrevMap)
library(lme4)

chikwawa <- read.csv("lab_data.csv")

attach(chikwawa) # attach variables to use

# Extraction of the elevation data
elev.data <- unique(chikwawa[,c("web_x","web_y","elevation")]) 

# Parameter estimation
fit.mle <- 
  linear.model.MLE(log(elevation) ~ 1, coords=~web_x+web_y,
                   data=elev.data,
                   start.cov.pars = 5.88,
                   fixed.rel.nugget = 0,
                   kappa=0.5) # this indicated the fitting of an exponential correlation

summary(fit.mle,log.cov.pars=FALSE)

variog.diagnostic.lm(fit.mle,param.uncertainty = TRUE,
                     which.test="variogram")

plot(elev.data[,c("web_x","web_y")],asp=1)

# Creation of the convex hull
ind.chull <- chull(elev.data[,c("web_x","web_y")]) 
poly <- elev.data[ind.chull,c("web_x","web_y")]
poly <- as.matrix(rbind(poly,poly[1,])) # closing the polygon
lines(poly)

# Creation of a 500 by 500 m prediction grid
library(splancs)
grid.pred <- gridpts(poly,xs=0.5,ys=0.5) 
points(grid.pred,pch=20,cex=0.1)

# Spatial prediction (log-elevation)
pred.mle <- 
spatial.pred.linear.MLE(fit.mle,grid.pred = grid.pred,
                        scale.predictions = "logit",
                        standard.errors = TRUE,
                        quantiles = c(0.025,0.975))

plot(pred.mle,type="logit",summary="predictions")
plot(pred.mle,type="logit",summary="standard.errors")
plot(pred.mle,type="logit",summary="quantiles",zlim=c(3.8,6.3))

# Spatial prediction (elevation)
pred.mle <- 
  spatial.pred.linear.MLE(fit.mle,grid.pred = grid.pred,
                          standard.errors = TRUE,
                          quantiles = c(0.025,0.975),
                          n.sim.prev=1000)

dim(pred.mle$samples)

elevation.samples <- exp(pred.mle$samples)

elevation.predictions <- apply(elevation.samples,1,mean)
elevation.se <- apply(elevation.samples,1,sd)
elevation.q025 <- apply(elevation.samples,1,function(x) quantile(x,0.025))
elevation.q975 <- apply(elevation.samples,1,function(x) quantile(x,0.975))

r.predicitons <- rasterFromXYZ(cbind(grid.pred,elevation.predictions))
r.se <- rasterFromXYZ(cbind(grid.pred,elevation.se))
r.q025 <- rasterFromXYZ(cbind(grid.pred,elevation.q025))
r.q975 <- rasterFromXYZ(cbind(grid.pred,elevation.q975))

plot(r.predicitons)
plot(r.se)

par(mfrow=c(1,2))
plot(r.q025,zlim=c(45,520))
plot(r.q975,zlim=c(45,520))
