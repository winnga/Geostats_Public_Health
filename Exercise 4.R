rm(list=ls()) 

library(PrevMap)
library(lme4)

chikwawa <- read.csv("lab_data.csv")

check.corr <- 
spat.corr.diagnostic(Pf ~ age.months,
                     coords=~web_x+web_y,
                     likelihood="Binomial",
                     ID.coords = chikwawa$ID,
                     data=chikwawa,
                     uvec=seq(0,20,length=15),
                     n.sim = 10000,
                     lse.variogram = TRUE)

glm.fit <- glm(Pf ~ age.months,data=chikwawa,family=binomial)

theta0 <- c(coef(glm.fit),check.corr$lse.variogram)

control.mcmc <- control.mcmc.MCML(n.sim=10000,burnin=2000,thin=8)

fit.MCML <- 
binomial.logistic.MCML(Pf ~ age.months,
                       units.m=~units.m,
                       coords=~web_x+web_y,
                       data = chikwawa,kappa=0.5,
                       par0=theta0,
                       ID.coords=chikwawa$ID,
                       start.cov.pars = c(52,0.5),
                       control.mcmc = control.mcmc)

summary(fit.MCML)

# Creation of the convex hull
coords <- unique(chikwawa[,c("web_x","web_y")])
ind.chull <- chull(coords) 
poly <- coords[ind.chull,]
poly <- as.matrix(rbind(poly,poly[1,])) # closing the polygon

# Creation of a 500 by 500 m prediction grid
library(splancs)
grid.pred <- gridpts(poly,xs=0.5,ys=0.5) 

predictors <- data.frame(age.months=rep(20,nrow(grid.pred)))
pred.Pf <- 
spatial.pred.binomial.MCML(fit.MCML,grid.pred = grid.pred,
                        predictors=predictors,
                        control.mcmc = control.mcmc,
                        scale.predictions = c("prevalence","odds"),
                        standard.errors=TRUE)

plot(pred.Pf,type="prevalence",summary="predictions")
plot(pred.Pf,type="prevalence",summary="standard.errors")

plot(pred.Pf,type="odds",summary="predictions")
plot(pred.Pf,type="odds",summary="standard.errors")
