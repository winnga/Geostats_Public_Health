rm(list=ls()) 

library(PrevMap)
source("test_version.R")
library(lme4)

chikwawa <- read.csv("lab_data.csv")

# Fitting of the linear model

fit.mle <- 
  linear.model.MLE(log(Hb) ~ sex, coords=~web_x+web_y,
                   kappa=0.5,ID.coords = chikwawa$ID,
                   start.cov.pars = c(4.8,0.02/0.038),data=chikwawa,
                   fixed.rel.nugget = 0)

# Creation of the convex hull
coords <- unique(chikwawa[,c("web_x","web_y")])
ind.chull <- chull(coords) 
poly <- coords[ind.chull,]
poly <- as.matrix(rbind(poly,poly[1,])) # closing the polygon

# Creation of a 500 by 500 m prediction grid
library(splancs)
grid.pred <- gridpts(poly,xs=0.5,ys=0.5) 


predictors <- data.frame(sex=factor(rep("Female",nrow(grid.pred)),levels = c("Female","Male")))
pred.Hb <- 
spatial.pred.linear.MLE(fit.mle,grid.pred = grid.pred,
                        predictors=predictors,
                        scale.predictions = "logit",
                        n.sim.prev = 1000)
plot()

omega2.hat <- coef(fit.mle)["omega^2"]
prev.samples <- pnorm((log(13.5)-pred.Hb$samples)/sqrt(omega2.hat))

prev.predictions <- apply(prev.samples,1,mean)
prev.se <- apply(prev.samples,1,sd)
prev.q025 <- apply(prev.samples,1,function(x) quantile(x,0.025))
prev.q975 <- apply(prev.samples,1,function(x) quantile(x,0.975))

r.predictions <- rasterFromXYZ(cbind(grid.pred,prev.predictions))
r.se <- rasterFromXYZ(cbind(grid.pred,prev.se))
r.q025 <- rasterFromXYZ(cbind(grid.pred,prev.q025))
r.q975 <- rasterFromXYZ(cbind(grid.pred,prev.q975))

plot(r.predictions)
plot(r.se)
