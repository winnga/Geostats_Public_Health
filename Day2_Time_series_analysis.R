rm(list=ls())

###################################################
###################################################
############### TEMPERATURE DATA ##################
###################################################
###################################################

lancs <- read.csv("maxtemp.csv")
lancs$t <- 1:nrow(lancs)
plot(lancs$t,lancs$temp,type="l",
     xlab="Day",ylab="Temperature (C)")

lm.fit <- lm(temp ~ sin(2*pi*t/365)+cos(2*pi*t/365),
             data=lancs)

lines(predict(lm.fit))

plot(lancs$t,residuals(lm.fit),type="l")

acf(residuals(lm.fit))

xreg.temp <- cbind(sin(2*pi*lancs$t/365),
                     cos(2*pi*lancs$t/365))
arima.fit <- 
arima(x=lancs$temp,xreg = xreg.temp,
      order=c(1,0,0))

arima.fit

acf(arima.fit$residuals)

temp.pred <- arima.fit$residuals+lancs$temp

matplot(lancs$t,cbind(lancs$temp,
                      temp.pred),type="l",
     xlab="Day",ylab="Temperature (C)",
     col=1:2)



####################################################
####################################################
################## PLAGUE DATA #####################
####################################################
####################################################

plague <- read.csv("plague.csv")
plague$t <- 1:nrow(plague)
str(plague)

plot(plague$t,plague$cases,type="l",
     xlab="Time",ylab="Cases")

plague$incidence <- plague$cases/plague$pop.at.risk

plot(plague$t,plague$incidence*100000,type="l",
     xlab="Time",ylab="Incidence (per 100,000)")

library(lme4)
glmer.fit <- glmer(cases ~ offset(log(pop.at.risk))+
                           sin(2*pi*t/12)+cos(2*pi*t/12)+
                           (1|t),
                   data=plague,family=poisson)
Z.hat <- ranef(glmer.fit)$t[,1]
plot(Z.hat,type="l")
acf(Z.hat)

residuals.data <- data.frame(Z.hat=Z.hat,
                             t=plague$t,
                             t2=1)

library(PrevMap)
plot(variogram(residuals.data,~Z.hat,~t+t2,
               uvec=seq(0,30,length=15)),
     type="l",
     xlab="Time separation (months)")

plague$t2 <- 1
spat.corr.diagnostic(cases ~ sin(2*pi*t/12)+cos(2*pi*t/12),
                     data=plague,
                     coords=~t+t2,
                     uvec=seq(0,30,length=15),
                     likelihood="Poisson",lse.variogram = TRUE)

mcml <- control.mcmc.MCML(n.sim=10000,burnin=2000,thin=8)

beta.guess <- glmer.fit@beta
sigma2.guess <- 0.31
phi.guess <- 9.32
tau2.guess <- 0.09

par0.plague <- c(beta.guess,sigma2.guess,phi.guess,tau2.guess)

fit.MCML <- 
poisson.log.MCML(cases~sin(2*pi*t/12)+cos(2*pi*t/12),
                 units.m=~pop.at.risk,
                 control.mcmc = mcml,
                 par0=par0.plague,kappa=0.5,
                 start.cov.pars = c(phi.guess,
                                    tau2.guess/sigma2.guess),
                 coords=~t+t2,
                 data=plague)

summary(fit.MCML)

variog.diagnostic.glgm(fit.MCML)

time.pred <- cbind(t=plague$t,t2=1)

pred.incidence <- 
spatial.pred.poisson.MCML(fit.MCML,
                          predictors = data.frame(t=plague$t),
                          grid.pred = time.pred,
                          scale.predictions = "exponential",
                          quantile=c(0.025,0.975),
                          control.mcmc = mcml)


incidence.hat <- pred.incidence$exponential$predictions
incidence.q025 <- pred.incidence$exponential$quantiles[,1]
incidence.q975 <- pred.incidence$exponential$quantiles[,2]

par(mfrow=c(1,1))
matplot(plague$t,
        cbind(plague$incidence,
              incidence.hat,
              incidence.q025,
              incidence.q975)*100000,
        col=c(1,2,2,2),
        lty=c("solid","solid","dashed","dashed"),
        type="l",
        xlab="Time",
        ylab="Incidence (per 100,000)")
