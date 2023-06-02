rm(list=ls())
library(PrevMap)

an <- read.csv("anopheles.csv")

str(an)

plot(an$elevation,log(an$An.gambiae+1),pch=20)
abline(glm(An.gambiae ~ elevation,data=an,family = poisson),
       col=2,lwd=2)


spat.corr.diagnostic(An.gambiae ~ elevation,
                     coords=~I(web_x/1000)+I(web_y/1000),
                     likelihood = "Poisson",
                     data=an,
                     uvec=seq(0,10,length=15),
                     n.sim = 10000)

library(lme4)
an$ID <- 1:nrow(an)
glmer.fit <- glmer(An.gambiae ~ I(elevation/1000)+(1|ID), 
                   family=poisson,data=an)
summary(glmer.fit)

an$Z <- ranef(glmer.fit)$ID[,1]
vari <- variogram(an,~Z,
                  ~I(web_x/1000)+I(web_y/1000),
                  uvec=seq(0,10,length=15))

eyefit(vari)

fit.LA <- glgm.LA(An.gambiae ~ I(elevation/1000),
                  coords=~I(web_x/1000)+I(web_y/1000),
                  kappa=0.5,
                  start.cov.pars = 3,
                  fixed.rel.nugget = 0,
                  family="Poisson",
                  data=an,
                  return.covariance = FALSE)

par0 <- coef(fit.LA)

control.mcmc <- control.mcmc.MCML(n.sim = 10000,
                                  burnin=2000,
                                  thin=8)

fit.MCML <- 
poisson.log.MCML(An.gambiae ~ I(elevation/1000),
                 coords=~I(web_x/1000)+I(web_y/1000),
                 kappa=0.5,
                 start.cov.pars = 1.5,
                 fixed.rel.nugget = 0,
                 par0=par0,control.mcmc = control.mcmc,
                 data=an,method="nlminb")

summary(fit.MCML)
variog.diagnostic.glgm(fit.MCML,uvec=seq(0,10,length=15),
                       n.sim = 1000)

elev <- read.csv("anopheles_elevation.csv")

plot(rasterFromXYZ(elev))
points(an[,c("web_x","web_y")],pch=20)

pred.MCML <- spatial.pred.poisson.MCML(fit.MCML,
                          predictors = elev,
                          grid.pred = elev[,c("web_x","web_y")],
                          control.mcmc = control.mcmc,
                          scale.predictions = "exponential")

plot(pred.MCML,type="exponential")

glm.fit <- glm(An.gambiae ~ I(elevation/1000), data = an,
               family = poisson)

pred.GLM <- predict(glm.fit,newdata=elev,type="response")

r.MCML <- rasterFromXYZ(cbind(elev[,c("web_x","web_y")],
                              pred.MCML$exponential$predictions))
r.GLM <- rasterFromXYZ(cbind(elev[,c("web_x","web_y")],pred.GLM))

par(mfrow=c(1,2))
plot(r.MCML,zlim=c(2,13.5))
plot(r.GLM,zlim=c(2,13.5))

summary(fit.MCML)
summary(glm.fit)
