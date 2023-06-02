rm(list=ls())
library(PrevMap)
library(INLA)
ec.ozone <- read.csv("eastcoast_ozone.csv")
ec.bndrs <- read.csv("eastcoast_bndrs.csv")
ec.grid <-  read.csv("eastcoast_grid.csv")

ec.ozone[,c("web_x","web_y")] <- ec.ozone[,c("web_x","web_y")]/1000
ec.grid[,c("web_x","web_y")] <- ec.grid[,c("web_x","web_y")]/1000
ec.bndrs <- ec.bndrs/1000
###################################################################################
max.vec <- function(v,val) sapply(v,function(x) max(0,x-val))
lm.fit <- lm(ozone ~ log(popdens+1)+I(max.vec(log(popdens+1),5)),
             data=ec.ozone)
beta.hat <- coef(lm.fit)
log.pd.set <- seq(0,10,length=1000)
log.pd.set.vars <- cbind(1,log.pd.set,max.vec(log.pd.set,5))
broken.sticks <- as.numeric(log.pd.set.vars%*%beta.hat)
std.errors <- sqrt(diag(log.pd.set.vars%*%vcov(lm.fit)%*%t(log.pd.set.vars)))
ci.95 <- cbind(broken.sticks-qnorm(0.975)*std.errors,
               broken.sticks+qnorm(0.975)*std.errors)

plot(log(ec.ozone$popdens+1),ec.ozone$ozone,
     xlab="log(Pop. density + 1)",
     ylab="Ozone concentration")
matlines(log.pd.set,cbind(ci.95,broken.sticks),
         lty=c("dashed","dashed","solid"),
         col=2,type="l",lwd=2)

library(GenKern)
kern <- KernSur(ec.ozone$web_x,ec.ozone$web_y,
                xgridsize = 200,
                ygridsize = 200)
r.dens <- rasterFromXYZ(cbind(expand.grid(kern$xords,kern$yords),as.numeric(kern$zden)))
log.dens <- log(extract(r.dens,ec.ozone[,c("web_x","web_y")]))
plot(log(ec.ozone$popdens+1),log.dens)

lm.fit <- lm(log.dens ~ log(popdens+1)+I(max.vec(log(popdens+1),5)),
             data=ec.ozone)
beta.hat <- coef(lm.fit)
log.pd.set <- seq(0,10,length=1000)
log.pd.set.vars <- cbind(1,log.pd.set,max.vec(log.pd.set,5))
broken.sticks <- as.numeric(log.pd.set.vars%*%beta.hat)
std.errors <- sqrt(diag(log.pd.set.vars%*%vcov(lm.fit)%*%t(log.pd.set.vars)))
ci.95 <- cbind(broken.sticks-qnorm(0.975)*std.errors,
               broken.sticks+qnorm(0.975)*std.errors)

plot(log(ec.ozone$popdens+1),log.dens,
     xlab="log(Pop. density + 1)",
     ylab="log(Sampling density)")
matlines(log.pd.set,cbind(ci.95,broken.sticks),
         lty=c("dashed","dashed","solid"),
         col=2,type="l",lwd=2)

lm.fit <- lm(ozone ~ log.dens,
             data=ec.ozone)
beta.hat <- coef(lm.fit)
log.pd.set <- seq(min(log.dens),max(log.dens),length=1000)
log.pd.set.vars <- cbind(1,log.pd.set)
broken.sticks <- as.numeric(log.pd.set.vars%*%beta.hat)
std.errors <- sqrt(diag(log.pd.set.vars%*%vcov(lm.fit)%*%t(log.pd.set.vars)))
ci.95 <- cbind(broken.sticks-qnorm(0.975)*std.errors,
               broken.sticks+qnorm(0.975)*std.errors)

plot(log.dens,ec.ozone$ozone,
     xlab="log(Sampling density)",
     ylab="Ozone concentration")
matlines(log.pd.set,cbind(ci.95,broken.sticks),
         lty=c("dashed","dashed","solid"),
         col=2,type="l",lwd=2)

mesh <- inla.mesh.2d(loc=ec.ozone[,c("web_x","web_y")],
                     max.edge=c(30,300),cutoff = 100,
                     boundary=inla.mesh.segment(ec.bndrs))

lm.fit <- lm(log.dens ~ log(popdens+1)+I(max.vec(log(popdens+1),5)),
             data=ec.ozone)
beta.sd <- coef(lm.fit)

lm.fit <- lm(ozone ~ log(popdens+1)+I(max.vec(log(popdens+1),5)),
             data=ec.ozone)
beta.lm <- coef(lm.fit)

par0.pd <- set.par.ps(q=3,p=3,intensity = c(beta.sd, 1.7917534, 320.2511335),
                      response = c(beta.lm,42.24544,237.27294,16.23562 ),
                      preferentiality.par =3.327936 )
control.mcmc <- control.mcmc.MCML(n.sim=10000,burnin = 2000,thin=8)
fit.ps.pd <- lm.ps.MCML(ozone~log(popdens+1)+I(max.vec(log(popdens+1),5)),
                        ~log(popdens+1)+I(max.vec(log(popdens+1),5)),
                        coords=~web_x+web_y,
                        data.response = ec.ozone,
                        data.intensity = ec.grid,
                        par0=par0.pd,control.mcmc=control.mcmc,
                        kappa1=1,kappa2=0.5,mesh=mesh,grid.intensity=ec.grid[,1:2],
                        method = "nlminb")
summary(fit.ps.pd)

par0.pd <- coef(fit.ps.pd)
pred.ps.pd <- spatial.pred.lm.ps(fit.ps.pd,
                                 predictors=data.frame(popdens=ec.grid$popdens),control.mcmc=control.mcmc)

par(mfrow=c(2,1))
# predicted sampling intensity
plot(pred.ps.pd,target=2,main="Sampling intensity")
lines(ec.bndrs)
# predicted ozone concentration
plot(pred.ps.pd,target=1,main="Ozone concentration")
lines(ec.bndrs)


