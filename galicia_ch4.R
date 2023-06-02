rm(list=ls())
library(PrevMap)

galicia <- read.csv("galicia.csv")
galicia.bndrs <- read.csv("galicia_bndrs.csv")

galicia <- galicia[galicia$survey==2000,]

point.map(galicia,~lead,~x+y,xlab="",ylab="")
lines(galicia.bndrs)

par(mfrow=c(1,2))

diag.cor <- spat.corr.diagnostic(log(lead)~1,coords=~I(x/1000)+I(y/1000),
                                 data=galicia,likelihood = "Gaussian")

fit.MLE <- linear.model.MLE(log(lead)~1,coords=~I(x/1000)+I(y/1000),
                            start.cov.pars = 28,data=galicia,
                            fixed.rel.nugget = 0,
                            kappa=0.5)


diag.fit <- variog.diagnostic.lm(fit.MLE)


par(mfrow=c(1,2),mar=c(2,2,2,2))
plot(diag.cor,main=c("(a)"),ylim=c(0.07,0.3),xlab="Distance (km)",ylab="")
plot(diag.fit,main=c("(b)"),ylim=c(0.07,0.3),xlab="Distance (km)",ylab="")


cp <- control.prior(beta.mean=0,beta.covar=10^10,
                    uniform.sigma2 = c(0,100),
                    uniform.phi = c(0,100))
c.mcmc <- control.mcmc.Bayes(n.sim=10000,burnin=2000,thin=8,
                             epsilon.S.lim = c(0.05,1),
                             L.S.lim = c(5,100),h.theta3=NULL,
                             start.nugget = NULL,linear.model = TRUE)
fit.Bayes <- linear.model.Bayes(log(lead)~1,coords=~I(x/1000)+I(y/1000),
                                data=galicia,control.mcmc = c.mcmc,
                                control.prior = cp,kappa=0.5)

trace.plot(fit.Bayes,"beta",1)
trace.plot(fit.Bayes,"sigma2",1)
dens.plot(fit.Bayes,"sigma2",1)


galicia.grid <- gridpts(as.matrix(galicia.bndrs)/1000,xs=5,ys=5)

pred.MLE <- spatial.pred.linear.MLE(fit.MLE,grid.pred = galicia.grid,
                                    scale.predictions="logit",
                                    standard.errors = TRUE)
pred.Bayes <- spatial.pred.linear.Bayes(fit.Bayes,grid.pred = galicia.grid,
                                        scale.predictions="logit",
                                        standard.errors = TRUE)
plot(pred.Bayes,type="logit",summary="standard.errors")
plot(pred.MLE,type="logit",summary="standard.errors")


par(mfrow=c(1,2))
plot(pred.Bayes$logit$predictions,pred.MLE$logit$predictions,
     xlab="Bayesian inference",ylab="Likelihood-based inference",
     main="Predictions")
abline(0,1,col=2,lwd=2)
plot(pred.Bayes$logit$standard.errors,pred.MLE$logit$standard.errors,
     xlab="Bayesian inference",ylab="Likelihood-based inference",
     main="Standard errors")
abline(0,1,col=2,lwd=2)
mean(pred.Bayes$logit$standard.errors > pred.MLE$logit$standard.errors)

geo.summary <- summary(fit.MLE)
tab.mle <-
  cbind(t(geo.summary$coefficients[,1:2]),
        geo.summary$coefficients[,1]-qnorm(0.975)*geo.summary$coefficients[,2],
        geo.summary$coefficients[,1]+qnorm(0.975)*geo.summary$coefficients[,2])
tab.mle <- rbind(tab.mle,
                 (cbind(
                   geo.summary$cov.pars[,1:2],
                   geo.summary$cov.pars[,1]-qnorm(0.975)*geo.summary$cov.pars[,2],
                   geo.summary$cov.pars[,1]+qnorm(0.975)*geo.summary$cov.pars[,2])))

fit.Bayes$estimate[,-1] <- log(fit.Bayes$estimate[,-1])
tab.Bayes <- cbind(apply(fit.Bayes$estimate,2,mean),
                   apply(fit.Bayes$estimate,2,sd),
                   apply(fit.Bayes$estimate,2,function(x) quantile(x,0.025)),
                   apply(fit.Bayes$estimate,2,function(x) quantile(x,0.975)))

tab <- rbind(tab.mle,tab.Bayes)

tab
