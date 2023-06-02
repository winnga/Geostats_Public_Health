rm(list = ls())

rb <- read.csv("LiberiaRemoData.csv")

beta.guess <- coef(glm(cbind(npos,ntest-npos) ~ log(elevation),
                  family=binomial,data=rb))

spat.corr.diagnostic(npos~log(elevation),
                     units.m=~ntest,data=rb,
                     coords=~I(utm_x/1000)+I(utm_y/1000),
                     likelihood = "Binomial",
                     lse.variogram = TRUE)

sigma2.guess <- 0.02
phi.guess <- 36.95
tau2.guess <-0.01

par0.rb <- c(beta.guess,sigma2.guess,phi.guess,tau2.guess)

mcml <- control.mcmc.MCML(n.sim=10000,burnin=2000,thin=8)

fit.mle <-
binomial.logistic.MCML(npos ~ log(elevation),
                       units.m = ~ ntest,
                       coords=~I(utm_x/1000)+I(utm_y/1000),
                       par0=par0.rb,control.mcmc = mcml,
                       kappa=0.5,
                       start.cov.pars = c(phi.guess,tau2.guess/sigma2.guess),
                       data=rb,method="nlminb")

summary(fit.mle)

par0.rb <- c(beta.guess,sigma2.guess,phi.guess)
fit.mle <-
  binomial.logistic.MCML(npos ~ I(log(elevation)),
                         units.m = ~ ntest,
                         coords=~I(utm_x/1000)+I(utm_y/1000),
                         par0=par0.rb,control.mcmc = mcml,
                         kappa=0.5,
                         fixed.rel.nugget = 0,
                         start.cov.pars = phi.guess,
                         data=rb,method="nlminb")


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
  
pred.mle <- 
spatial.pred.binomial.MCML(fit.mle,grid.pred = liberia.grid,
                           predictors = predictors.rb,
                           control.mcmc = mcml,
                           scale.predictions = c(
                             "logit","prevalence"),
                           thresholds = 0.2,
                           scale.thresholds = "prevalence")

plot(pred.mle,"prevalence","predictions")
plot(pred.mle,summary="exceedance.prob")

predictors.rb0 <- data.frame(elevation=rep(1,nrow(liberia.grid)))

pred.S <- 
  spatial.pred.binomial.MCML(fit.mle,grid.pred = liberia.grid,
                             predictors = predictors.rb0,
                             control.mcmc = mcml,
                             scale.predictions = "logit")

elev.contr <- rasterFromXYZ(cbind(liberia.grid,
                            pred.mle$logit$predictions-
                            pred.S$logit$predictions+
                              coef(fit.mle)[1]))


plot(pred.mle,"logit","predictions",zlim=c(-3,-0.5),
     main="Log(elevation)+S")
plot(elev.contr,zlim=c(-3,-0.5),main="Log(elevation)")

plot(elev.contr)

variog.diagnostic.glgm(fit.mle,
                       uvec=seq(20,150,length=15),
                       n.sim = 1000)
