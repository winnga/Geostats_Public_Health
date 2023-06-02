rm(list=ls())

library(PrevMap)
library(lme4)

chikwawa <- read.csv("lab_data.csv")
chikwawa$Hb.below12 <- 1*(chikwawa$Hb<12)

check.corr <- spat.corr.diagnostic(Hb.below12 ~ 1, coords = ~web_x+web_y,
                                   data=chikwawa,ID.coords = chikwawa$ID,
                                   likelihood = "Binomial",n.sim=10000,
                                   uvec=seq(2,12,length=15),
                                   lse.variogram = TRUE)
beta0 <- coef(glm(IRS ~ 1, data=chikwawa,family="binomial"))
cov.pars0 <- check.corr$lse.variogram[c("sigma^2","phi")]

theta0 <- c(beta0,cov.pars0)

control.mcmc <- control.mcmc.MCML(n.sim=10000,burnin=2000,thin=8)


fit.MCML <- 
binomial.logistic.MCML(Hb.below12 ~ 1, units.m = ~units.m,
                       coords=~web_x+web_y,
                       control.mcmc=control.mcmc,
                       par0=theta0,start.cov.pars = 1.8,
                       fixed.rel.nugget = 0,
                       ID.coords=chikwawa$ID,data=chikwawa,
                       kappa=0.5)

theta0 <- coef(fit.MCML)

fit.MCML <- 
  binomial.logistic.MCML(Hb.below12 ~ 1, units.m = ~units.m,
                         coords=~web_x+web_y,
                         control.mcmc=control.mcmc,
                         par0=theta0,start.cov.pars = 1.8,
                         fixed.rel.nugget = 0,
                         ID.coords=chikwawa$ID,data=chikwawa,
                         kappa=0.5)
summary(fit.MCML,log.cov.pars=FALSE)



# Creation of the convex hull
coords <- unique(chikwawa[,c("web_x","web_y")])
ind.chull <- chull(coords) 
poly <- coords[ind.chull,]
poly <- as.matrix(rbind(poly,poly[1,])) # closing the polygon
plot(poly,type="l")

# Creation of a 500 by 500 m prediction grid
library(splancs)
grid.pred <- gridpts(poly,xs=0.5,ys=0.5) 
points(grid.pred,pch=20,cex=0.1)

pred.MCML <- 
spatial.pred.binomial.MCML(fit.MCML,grid.pred=grid.pred,
                           control.mcmc = control.mcmc,
                           scale.predictions = c("prevalence","odds"),
                           standard.errors = TRUE)

plot(pred.MCML,type="prevalence",summary="predictions")
plot(pred.MCML,type="prevalence",summary="standard.errors")

plot(pred.MCML,type="odds",summary="predictions")
plot(pred.MCML,type="odds",summary="standard.errors")
