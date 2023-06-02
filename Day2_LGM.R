rm(list = ls())
library(PrevMap)

galicia <- read.csv("galicia.csv")
galicia <- galicia[galicia$survey==2000,]
galicia.bndrs <- read.csv("galicia_bndrs.csv")

point.map(galicia,~log(lead),coords=~x+y,
          pt.divide="quintiles")
lines(galicia.bndrs,type="l")

vari <- variogram(galicia,~lead,
          coords=~I(x/1000)+I(y/1000),
          uvec=seq(10,150,length=15))
plot(vari,type="b")

eyefit(vari)

spat.corr.diagnostic(lead~1, 
                     data=galicia,
                     coords=~I(x/1000)+I(y/1000),
                     likelihood = "Gaussian",
                     ID.coords = 1:nrow(galicia),
                     lse.variogram = TRUE,
                     uvec=seq(10,150,length=15))

sigma2.start <- 1.513003e+00
phi.start <- 1.710739e+01

lgm.fit.mle <- 
linear.model.MLE(lead~1, coords=~I(x/1000)+I(y/1000), 
                 data = galicia,
                 start.cov.pars = phi.start,
                 fixed.rel.nugget = 0,
                 kappa=0.5,method = "nlminb")

summary(lgm.fit.mle)

library(splancs)
grid.pred.galicia <- gridpts(as.matrix(galicia.bndrs),
                             xs=5000,ys=5000)/1000

pred.lead.MLE <- spatial.pred.linear.MLE(lgm.fit.mle,
                                     grid.pred = grid.pred.galicia,
                                     standard.errors = TRUE,
                                     scale.predictions = "logit",
                                     n.sim.prev = 1000)

plot(pred.lead.MLE,type="logit",summary="predictions")
plot(pred.lead.MLE,type="logit",summary="standard.errors")

max.samples.MLE <- apply(pred.lead.MLE$samples,2,max)
hist(max.samples.MLE)

####

log.prior.sigma2.galicia <- function(sigma2) {
  dunif(sigma2,0,10,log=TRUE)
}

log.prior.phi.galicia <- function(phi) {
  dunif(phi,0,50,log=TRUE)
}

cp <- 
control.prior(beta.mean = 0,beta.covar = 10^10,
              log.prior.sigma2 = log.prior.sigma2.galicia,
              log.prior.phi = log.prior.phi.galicia)

bayes.mcmc <- 
control.mcmc.Bayes(n.sim=2000,burnin=1000,thin=1,
                   epsilon.S.lim = 0.1,
                   start.nugget = NULL,
                   L.S.lim = 5)

lgm.fit.bayes <- 
linear.model.Bayes(lead~1, coords=~I(x/1000)+I(y/1000),
                   data=galicia,kappa=0.5,
                   control.prior = cp,
                   control.mcmc = bayes.mcmc)

trace.plot(lgm.fit.bayes,"beta",component.beta = 1)
autocor.plot(lgm.fit.bayes,"beta",component.beta = 1)

trace.plot(lgm.fit.bayes,"sigma2")
autocor.plot(lgm.fit.bayes,"sigma2",component.beta = 1)

trace.plot(lgm.fit.bayes,"phi")
autocor.plot(lgm.fit.bayes,"phi",component.beta = 1)

pred.lead.Bayes <- 
spatial.pred.linear.Bayes(lgm.fit.bayes,
                            grid.pred = grid.pred.galicia,
                            scale.predictions = "logit",
                          standard.errors = TRUE)
plot(pred.lead.Bayes,"logit","predictions")
plot(pred.lead.Bayes,"logit","standard.errors")

max.samples.Bayes <- apply(pred.lead.Bayes$samples,1,max)
hist(max.samples.Bayes)

###

plot(pred.lead.MLE$logit$predictions,
     pred.lead.Bayes$logit$predictions,
     xlab="MLE",ylab="Bayes",
     main="Predictions")
abline(0,1,col=2,lwd=2)

plot(pred.lead.MLE$logit$standard.errors,
     pred.lead.Bayes$logit$standard.errors,
     xlab="MLE",ylab="Bayes",
     main="Predictions")
abline(0,1,col=2,lwd=2)

plot(density(max.samples.MLE),main="")
lines(density(max.samples.Bayes),col=2)

###

library(help=PrevMap)
variog.diagnostic.lm(lgm.fit.mle)
