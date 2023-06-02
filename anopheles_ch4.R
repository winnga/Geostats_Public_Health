rm(list=ls())
library(PrevMap)
mosq <- read.csv("anopheles.csv")

fit.LA <- glgm.LA(An.gambiae~elevation,
                  coords=~I(web_x/1000)+I(web_y/1000),kappa=0.5,
                  start.cov.pars = 20,fixed.rel.nugget = 0,
                  data=mosq,family="Poisson")

par0 <- coef(fit.LA)
c.mcmc <- control.mcmc.MCML(n.sim=42000,burnin=2000,thin=8)
fit.MCML <- poisson.log.MCML(An.gambiae~elevation,control.mcmc = c.mcmc,
                             par0=par0,
                             coords=~I(web_x/1000)+I(web_y/1000),kappa=0.5,
                             start.cov.pars = 1.523717640 ,fixed.rel.nugget = 0,
                             data=mosq,method="nlminb")

summary(fit.MCML)
