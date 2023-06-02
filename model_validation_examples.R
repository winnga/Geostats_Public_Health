rm(list=ls())
library(PrevMap)
coords <- expand.grid(seq(0,1,length=20),seq(0,1,length=20))
U <- as.matrix(dist(coords))

### Scenario 1: omission of the nugget effect
set.seed(1234)
sigma2 <- 1
phi <- 0.2
tau2 <- 1

Sigma <- sigma2*exp(-U/phi)
diag(Sigma) <- diag(Sigma)+tau2
mu <- 1
n <- nrow(coords)
y <- mu+t(chol(Sigma))%*%rnorm(n)
data.sim <- data.frame(x1=coords[,1],x2=coords[,2],y=y)

fit <- linear.model.MLE(y~1,data=data.sim,coords=~x1+x2,
                        fixed.rel.nugget = 0,kappa = 0.5,
                        start.cov.pars = 0.1,method = "nlminb")
summary(fit)
fit.diag <- variog.diagnostic.lm(fit,range.fact = 1,
                                 uvec=seq(0,0.5,length=15))


fit.correct <- linear.model.MLE(y~1,data=data.sim,coords=~x1+x2,
                                kappa = 0.5,
                                start.cov.pars = c(0.1,1),method = "nlminb")
fit.diag.correct <- variog.diagnostic.lm(fit.correct,range.fact=1,
                                         uvec=seq(0,0.5,length=15))

### Scenario 2: Miss-specification of the smoothness parameter
phi <- 0.01
tau2 <- 0

n <- nrow(coords)

sigma2 <- 1
Sigma <- matern(U,phi,5)
diag(Sigma) <- diag(Sigma)+tau2
mu <- 1

y <- mu+t(chol(Sigma))%*%rnorm(n)
data.sim <- data.frame(x1=coords[,1],x2=coords[,2],y=y)

fit <- linear.model.MLE(y~1,data=data.sim,coords=~x1+x2,
                        kappa = 0.5,fixed.rel.nugget = 0,
                        start.cov.pars = 0.1,method = "nlminb")
fit.diag <- variog.diagnostic.lm(fit,range.fact=1)



fit.correct <- linear.model.MLE(y~1,data=data.sim,coords=~x1+x2,
                                kappa = 5,fixed.rel.nugget = 0,
                                start.cov.pars = 0.01,method = "nlminb")
fit.diag.correct <- variog.diagnostic.lm(fit.correct,range.fact = 1)


### Scenario 3: non-Gaussian data
set.seed(1234)
phi <- 0.2
tau2 <- 1

n <- nrow(coords)

sigma2 <- 2
Sigma <- sigma2*matern(U,phi,0.5)
diag(Sigma) <- diag(Sigma)+tau2
mu <- 1

y <- exp(mu+t(chol(Sigma))%*%rnorm(n))
data.sim <- data.frame(x1=coords[,1],x2=coords[,2],y=y)

fit <- linear.model.MLE(y~1,data=data.sim,coords=~x1+x2,
                        kappa = 0.5,
                        start.cov.pars = c(0.2,0.5),method = "nlminb")

fit.diag <- variog.diagnostic.lm(fit,range.fact=1)

fit.correct <- linear.model.MLE(log(y)~1,data=data.sim,coords=~x1+x2,
                                kappa = 0.5,
                                start.cov.pars = c(0.2,0.5),method = "nlminb")

fit.correct <- variog.diagnostic.lm(fit.correct,range.fact=1)



