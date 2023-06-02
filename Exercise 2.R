rm(list=ls()) 

library(PrevMap)
library(lme4)

chikwawa <- read.csv("lab_data.csv")

str(chikwawa)

attach(chikwawa) # attach variables to use

# Getting the initial values
check.corr <- 
spat.corr.diagnostic(log(Hb) ~ sex, coords=~web_x+web_y,
                     kappa=0.5,ID.coords = chikwawa$ID,
                     data=chikwawa,likelihood = "Gaussian",
                     which.test = "variogram",lse.variogram = TRUE)

fit.lmer <- 
lmer(log(Hb) ~ sex + (1|ID),data=chikwawa)

Z.hat <- ranef(fit.lmer)$ID[,1]
omega2.start <- var(log(chikwawa$Hb)-Z.hat[chikwawa$ID])
omega2.start

# Fitting of the linear model

fit.mle <- 
linear.model.MLE(log(Hb) ~ sex, coords=~web_x+web_y,
                 kappa=0.5,ID.coords = chikwawa$ID,
                 start.cov.pars = c(4.8,0.02/0.038),data=chikwawa,
                 fixed.rel.nugget = 0)

summary(fit.mle,log.cov.pars=FALSE)                 

# Model validation
diag.lm <- variog.diagnostic.lm(fit.mle,which.test = "variogram",uvec=seq(1,25,1))
sigma2.hat <- coef(fit.mle)["sigma^2"]
phi <- coef(fit.mle)["phi"]
curve(sigma2.hat*(1-exp(-x/phi.hat)),add=TRUE,col=2)
