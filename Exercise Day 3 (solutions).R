rm(list = ls())

kericho <- read.csv("Kericho.csv")

kericho$t <- 1:nrow(kericho)

plot(kericho$t,log(kericho$Cases),type="l",
     xlab="Month",ylab="log(Cases)")

max.vec <- function(vec,c) sapply(vec,function(x) max(x,c))

lm.fit <- lm(log(Cases) ~ I(t>230)+t+max.vec(t,60),
             data=kericho)
lines(predict(lm.fit),col=2)
summary(lm.fit)

plot(residuals(lm.fit),type="l")
acf(residuals(lm.fit),lag.max = 100)

lm.fit2 <- lm(log(Cases) ~ I(t>230)+t+max.vec(t,60)+
                           sin(2*pi*t/12)+cos(2*pi*t/12),
             data=kericho)

plot(residuals(lm.fit2),type="l")
acf(residuals(lm.fit2),lag.max = 100)


lm.fit3 <- lm(log(Cases) ~ I(t>230)+t+max.vec(t,60)+
                sin(2*pi*t/12)+cos(2*pi*t/12)+
                sin(2*pi*t/6)+cos(2*pi*t/6),
              data=kericho,x=TRUE)
plot(residuals(lm.fit3),type="l")
acf(residuals(lm.fit3),lag.max = 100)

plot(kericho$t,kericho$Cases,type="l",
     xlab="Month",ylab="log(Cases)")
lines(exp(predict(lm.fit3)),col=2)


xreg.kericho <- lm.fit3$x[,-1]

arima.fit <- 
arima(x=log(kericho$Cases), 
      xreg = xreg.kericho,
      order=c(1,0,0))

arima.fit

log.cases.pred <- arima.fit$residuals+log(kericho$Cases)

matplot(kericho$t,
     cbind(log(kericho$Cases),
             log.cases.pred),type="l",
     xlab="Month",ylab="log(Cases)",
     lty="solid")

acf(arima.fit$residuals)


#########
library(lme4)
glmer.fit <- glmer(Cases ~ I(t>230)+t+max.vec(t,60)+
                     sin(2*pi*t/12)+cos(2*pi*t/12)+
                     sin(2*pi*t/6)+cos(2*pi*t/6)+(1|t),
                   data=kericho,family=poisson) 
kericho$t2 <- 1
spat.corr.diagnostic(Cases ~ I(t>230)+t+max.vec(t,60)+
                       sin(2*pi*t/12)+cos(2*pi*t/12)+
                       sin(2*pi*t/6)+cos(2*pi*t/6),
                     coords=~t+t2,likelihood = "Poisson",
                     uvec=seq(0,20,length=15),
                     data=kericho,
                     lse.variogram = TRUE)

library(PrevMap)
mcml <- control.mcmc.MCML(n.sim=10000,burnin=2000,thin=8)

beta.guess <- glmer.fit@beta
sigma2.guess <- 2
phi.guess <- 1.11

par0.kericho <- c(beta.guess,sigma2.guess,phi.guess)

fit.MCML <- poisson.log.MCML(Cases ~ I(t>230)+t+max.vec(t,60)+
                               sin(2*pi*t/12)+cos(2*pi*t/12)+
                               sin(2*pi*t/6)+cos(2*pi*t/6),
                             data=kericho,
                             coords=~t+t2,
                             kappa=0.5,
                             fixed.rel.nugget = 0,
                             control.mcmc = mcml,
                             par0=par0.kericho,
                             start.cov.pars = phi.guess,
                             method="nlminb")


time.pred <- data.frame(t=kericho$t+0.01,t2=1)

pred.cases <- 
spatial.pred.poisson.MCML(fit.MCML,
                          grid.pred = time.pred,
                          control.mcmc = mcml,
                          predictors = data.frame(t=kericho$t),
                          scale.predictions = "exponential",
                          quantiles = c(0.025,0.975))

matplot(kericho$t,
        cbind(kericho$Cases,
              pred.cases$exponential$predictions,
              pred.cases$exponential$quantiles[,1],
              pred.cases$exponential$quantiles[,2]),
        col=c(1,2,2,2),
        lty=c("solid","solid","dashed","dashed"),
        type="l",
        xlab="Month",ylab="Cases")
