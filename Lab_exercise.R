rm(list=ls())

rb <- read.csv(file.choose())

str(rb)
rb$logit <- log((rb$npos+0.5)/(rb$ntest+0.5-rb$npos))

hist(rb$logit)

spat.corr.diagnostic(logit ~ 1,
                     data=rb, likelihood = "Gaussian",
                     coords=~I(utm_x/1000)+
                             I(utm_y/1000),
                     uvec=seq(0,200,length=15))


lm.fit <- lm(logit ~ 1, data=rb)
rb$residuals <- residuals(lm.fit)

vari <- variogram(rb,var.name = ~residuals,
                  coords = ~I(utm_x/1000)+
                            I(utm_y/1000),
                  uvec=seq(0,200,length=15))
eyefit(vari)

#sigma2=1.122
#phi=160
#tau2=0.14

fit.mle <- 
linear.model.MLE(logit ~ 1, 
                 start.cov.pars = c(160,0.14/1.122),
                 coords = ~I(utm_x/1000)+
                   I(utm_y/1000),
                 kappa=0.5,
                 data=rb)


summary(fit.mle)


liberia <- read.csv(file.choose()) # Liberia_bndrs.csv
liberia <- liberia[,-1]
liberia <- liberia/1000

plot(liberia,type="l")


library(splancs)
grid <- gridpts(as.matrix(liberia),xs=10,ys=10)

points(grid,pch=20)

pred <- 
spatial.pred.linear.MLE(fit.mle,grid.pred = grid,
                        scale.predictions = "prevalence",
                        n.sim.prev = 1000)

plot(pred,"prevalence","predictions")
