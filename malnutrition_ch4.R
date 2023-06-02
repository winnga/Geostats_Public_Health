rm(list=ls())

library(PrevMap)
Ghana.bndrs <- read.csv("Ghana_bndrs.csv")

maln <- read.csv("malnutrition.csv")

ID.coords <- create.ID.coords(maln,~utm_x+utm_y)

max.vec <- function(const,x) sapply(x, function(x.fun) max(const,x.fun))

maln.diag <- spat.corr.diagnostic(formula=HAZ ~ age+I(max.vec(1,age))+I(max.vec(2,age))+edu+wealth,
                                  coords=~utm_x+utm_y,ID.coords = ID.coords,
                                  data=maln,likelihood="Gaussian",
                                  uvec=seq(0,600,length=15),
                                  lse.variogram = TRUE)

geo.fit <- linear.model.MLE(HAZ ~age+I(max.vec(1,age))+I(max.vec(2,age))+
                              edu+wealth,
                            coords=~utm_x+utm_y, ID.coords=ID.coords,
                            start.cov.pars=c(20,2),kappa=0.5,data=maln,
                            fixed.rel.nugget = 0,
                            method="nlminb")
beta.hat.geo <- summary(geo.fit)$coefficients[,1:2]

cbind(beta.hat.geo,
             beta.hat.geo[,1]-qnorm(0.975)*beta.hat.geo[,2],
             beta.hat.geo[,1]+qnorm(0.975)*beta.hat.geo[,2])


library(splancs)
Ghana.grid <- gridpts(as.matrix(Ghana.bndrs),xs=10,ys=10)

Ghana.predictors <- data.frame(age=rep(2,nrow(Ghana.grid)),
                               edu=rep(1,nrow(Ghana.grid)),
                               wealth=rep(1,nrow(Ghana.grid)))

geo.pred <- spatial.pred.linear.MLE(geo.fit,grid.pred = Ghana.grid,
                                    predictors = Ghana.predictors,
                                    scale.predictions = "logit",
                                    standard.errors = TRUE,
                                    thresholds = -2,
                                    scale.thresholds = "logit")
geo.pred$exceedance.prob <- 1-geo.pred$exceedance.prob

par(mfrow=c(3,1),mar=c(2,2,3,2))
plot(geo.pred,type="logit",summary="predictions",main="(a) HAZ estimates")
lines(Ghana.bndrs)
plot(geo.pred,type="logit",summary="standard.errors",main="(b) Standard errors")
lines(Ghana.bndrs)
plot(geo.pred,summary="exceedance.prob",zlim=c(0,1),main="(c) Stunting risk")
lines(Ghana.bndrs)

geo.fit.diag <- variog.diagnostic.lm(geo.fit)

sample.age <- function(n.sim) {
  t.age <- log(maln$age/5/(1-maln$age/5))
  bw <- density(t.age)$bw
  t.age.obs <- t.age[sample(1:nrow(maln),n.sim,replace=TRUE)]
  t.age.sim <- rnorm(n.sim,mean=t.age.obs,sd=bw)
  5*exp(t.age.sim)/(1+exp(t.age.sim))
}

sample.edu.wealth <- function(n.sim) {
  tab <- table(maln$edu,maln$wealth)/nrow(maln)
  values <- expand.grid(unique(maln$edu),unique(maln$wealth))
  out <- values[sample(1:nrow(values),n.sim,prob=as.numeric(tab),replace = TRUE),]
  data.frame(edu=out[,1],wealth=out[,2])
}

n.sim <- 1000

predictors.samples <- list()
age.samples <- sample.age(n.sim)
edu.wealth.samples <- sample.edu.wealth(n.sim)
for(i in 1:n.sim) {
  predictors.samples[[i]] <- data.frame(age=rep(age.samples[i],nrow(Ghana.grid)),
                                        edu=rep(edu.wealth.samples[i,1],nrow(Ghana.grid)),
                                        wealth=rep(edu.wealth.samples[i,2],nrow(Ghana.grid)))
}

geo.pred.pop <- spatial.pred.linear.MLE(geo.fit,grid.pred = Ghana.grid,
                                    predictors.samples = predictors.samples,
                                    scale.predictions = "logit",
                                    type="joint",
                                    standard.errors = TRUE,
                                    thresholds = -2,n.sim.prev = n.sim,
                                    scale.thresholds = "logit")
geo.pred.pop$exceedance.prob <- 1-geo.pred$exceedance.prob

par(mfrow=c(3,1),mar=c(2,2,3,2))
plot(geo.pred.pop,type="logit",summary="predictions",main="(a) HAZ estimates")
lines(Ghana.bndrs)
plot(geo.pred.pop,type="logit",summary="standard.errors",main="(b) Standard errors")
lines(Ghana.bndrs)
plot(geo.pred.pop,summary="exceedance.prob",zlim=c(0,1),main="(c) Stunting risk")
lines(Ghana.bndrs)
