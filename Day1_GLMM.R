rm(list=ls())

library(lme4)

rb <- read.csv("LiberiaRemoData.csv")

str(rb)

rb$prev <- rb$npos/rb$ntest

plot(rb$elevation,rb$prev,
     xlab="Elevation", ylab="Prevalence")
plot(log(rb$elevation),rb$prev,
     xlab="log(Elevation)", ylab="Prevalence")

rb$e.logit <- log((rb$npos+0.5)/(rb$ntest-rb$npos+0.5))
plot(log(rb$elevation),rb$e.logit,
     xlab="log(Elevation)",
     ylab="Empirical logit")

max.vec <- function(vec,c) sapply(vec,function(x) max(x,c))
glm.fit <- glm(cbind(npos,ntest-npos) ~ log(elevation) + 
               max.vec(log(elevation),5.5),
               family=binomial, data=rb)

summary(glm.fit)

elevation.values <- seq(min(rb$elevation),
                        max(rb$elevation),
                        length=1000)
pred.val <- predict(glm.fit,newdata = data.frame(elevation=elevation.values))

lines(log(elevation.values),pred.val,col=2)

n.sim <- 1000
glm.fit0 <- glm(cbind(npos,ntest-npos) ~ log(elevation),
                family=binomial, data=rb, x=TRUE)
pred.val0 <- predict(glm.fit0,newdata = data.frame(elevation=elevation.values))
lines(log(elevation.values),pred.val0,col=3)

D.obs <- 2*(logLik(glm.fit)-logLik(glm.fit0))

1-pchisq(D.obs,1)

prevalence.val0 <- predict(glm.fit0,
                   newdata = data.frame(elevation=elevation.values),
                   type="response")
plot(rb$elevation,rb$prev,
     xlab="Elevation", ylab="Prevalence")
lines(elevation.values,prevalence.val0,col=2)


Sigma.par <- vcov(glm.fit0)
Sigma.par.sroot <- t(chol(Sigma.par))
beta.hat <- coef(glm.fit0) 
p <- length(beta.hat)
sim <- 
sapply(1:n.sim,function(i) {
                 beta.sample <- beta.hat+Sigma.par.sroot%*%rnorm(p)
                 lp.sample <- beta.sample[1]+beta.sample[2]*log(elevation.values)
                 prev.sample <- exp(lp.sample)/(1+exp(lp.sample))
               })
CI <- apply(sim,1,function(x) quantile(x,c(0.025,0.975)))

matlines(elevation.values,
         t(CI),col=2,type="l",lty="dashed")

#####
rb$ID <- 1:nrow(rb)
glmer.fit <- glmer(cbind(npos,ntest-npos) ~ log(elevation)+(1|ID),
                   data=rb,family=binomial,
                   nAGQ=50)

summary(glmer.fit)
summary(glm.fit0)

beta.hat.glmer <- glmer.fit@beta
lp.val.glmer <- beta.hat.glmer[1]+beta.hat.glmer[2]*log(elevation.values)
prevalence.val.glmer <- exp(lp.val.glmer)/(1+exp(lp.val.glmer))
vcov(glmer.fit)

Sigma.glmer.par <- vcov(glmer.fit)
Sigma.glmer.par.sroot <- t(chol(Sigma.glmer.par))

sim <- 
  sapply(1:n.sim,function(i) {
    beta.sample <- beta.hat.glmer+Sigma.glmer.par.sroot%*%rnorm(p)
    lp.sample <- beta.sample[1]+beta.sample[2]*log(elevation.values)
    prev.sample <- exp(lp.sample)/(1+exp(lp.sample))
  })
CI.glmer <- apply(sim,1,function(x) quantile(x,c(0.025,0.975)))

lines(elevation.values,prevalence.val.glmer,col=3)
matlines(elevation.values,
         t(CI.glmer),col=3,type="l",lty="dashed")
