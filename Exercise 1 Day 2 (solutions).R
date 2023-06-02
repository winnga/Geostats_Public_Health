rm(list = ls())

rb <- read.csv("LiberiaRemoData.csv")
rb$logit <- log((rb$npos+0.5)/(rb$ntest-rb$npos+0.5))


spat.corr.diagnostic(logit~log(elevation),
                     data=rb,
                     coords=~I(utm_x/1000)+I(utm_y/1000),
                     likelihood = "Gaussian",
                     lse.variogram = TRUE)

sigma2.guess <- 0.26
phi.guess <- 75.83
tau2.guess <-0.19


fit.mle <- linear.model.MLE(logit ~ log(elevation),
                            coords=~I(utm_x/1000)+I(utm_y/1000),
                            kappa=0.5,
                            start.cov.pars = c(phi.guess,tau2.guess/sigma2.guess),
                         data=rb,method="nlminb")

summary(fit.mle)

liberia.bndrs <- read.csv("Liberia_bndrs.csv")
library(splancs)  
liberia.grid <- gridpts(as.matrix(liberia.bndrs[,c("utm_x","utm_y")])/1000,
                        xs=3,ys=3)  

elevation <- raster("./LBR_alt/LBR_alt.gri")
elevation <- projectRaster(elevation,crs=CRS("+init=epsg:32629"))
elevation.pred <- extract(elevation,liberia.grid*1000)

ind.na <- which(is.na(elevation.pred))
liberia.grid <- liberia.grid[-ind.na,]
elevation.pred <- elevation.pred[-ind.na]

predictors.rb <- data.frame(elevation=elevation.pred)

pred.mle.lm <- 
  spatial.pred.linear.MLE(fit.mle,grid.pred = liberia.grid,
                          predictors = predictors.rb,
                          scale.predictions = "prevalence",
                          thresholds = 0.2,n.sim.prev = 1000,
                          scale.thresholds = "prevalence")

plot(pred.mle.lm,"prevalence","predictions")
plot(pred.mle.lm,summary="exceedance.prob")

load("binomial_predictions_rb.RData")

plot(pred.mle,"prevalence","predictions")
plot(pred.mle,summary="exceedance.prob")

plot(pred.mle.lm$prevalence$predictions,
     pred.mle$prevalence$predictions,
     xlab="Linear model (empirical logit)",
     ylab="Binomial model",pch=20)
abline(0,1,col=2,lwd=2)

plot(pred.mle.lm$exceedance.prob,
     pred.mle$exceedance.prob,
     xlab="Linear model (empirical logit)",
     ylab="Binomial model",pch=20)
abline(0,1,col=2,lwd=2)

r.diff.pred <- rasterFromXYZ(cbind(liberia.grid,
               pred.mle$prevalence$predictions-
               pred.mle.lm$prevalence$predictions)) 
plot(r.diff.pred)


r.diff.ep <- rasterFromXYZ(cbind(liberia.grid,
                                pred.mle$exceedance.prob-
                                pred.mle.lm$exceedance.prob)) 
plot(r.diff.ep)
