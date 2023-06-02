rm(list=ls())
library(PrevMap)

# To install INLA run the following line of code
# install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/testing")
library(INLA)


# Load the data
data("galicia")

# Load the boundary of Galicia
data("galicia.boundary")

# Re-scale geographical coordinates so that
# distances are in km
galicia$x <- galicia$x/1e3
galicia$y <- galicia$y/1e3
galicia.boundary <- galicia.boundary/1e3

# Geographical coordinates for the two surveys
coords97 <- galicia[galicia$survey==1997,c("x","y")]
coords00 <- galicia[galicia$survey==2000,c("x","y")]

# Quintiles plot of lead pollution in the two surveys
lead.quintiles <- quantile(galicia$lead,seq(0,1,0.2))
lead.class <- as.numeric(cut(galicia$lead,
                             breaks=lead.quintiles,include.lowest=TRUE))
par(mfrow=c(1,2))
plot(galicia.boundary,type="l",asp=1,xlab="",ylab="" )
points(coords97,pch=lead.class[galicia$survey==1997])
plot(galicia.boundary,type="l",asp=1,xlab="",ylab="" )
points(coords00,pch=lead.class[galicia$survey==2000])
legend(500,4980,c("1st quint.","2nd quint.", "3rd quint.","4th quint.", "5th quint."),
       pch=1:5,cex=1)

library(GenKern)
kern <- KernSur(coords97[,1],coords97[,2],
                range.y=c(4628.86, 4849.553),
                xgridsize = 200,
                ygridsize = 200)
r.dens <- rasterFromXYZ(cbind(expand.grid(kern$xords,kern$yords),as.numeric(kern$zden)))
poly <- SpatialPolygons(list(Polygons(list(Polygon(galicia.boundary)),1)))
r.dens <- mask(r.dens,poly)

par(mfrow=c(1,2),mar=c(4,4,2,6))
plot(r.dens)     
lines(poly)
log.dens <- log(extract(r.dens,coords97))
plot(log.dens,log(galicia$lead[galicia$survey==1997]),
     pch=20,
     xlab="log(sampling density)",
     ylab="log(Lead concentration)")
abline(lm(log(galicia$lead[galicia$survey==1997]) ~ log.dens))


# Set the parameter of the importance sampling distribution
par0 <- set.par.ps(p=2,q=1,intensity = c(alpha=-0.356,sigma2_1=6.610, phi1=185.5),
                   response = c(mu97=0.854, mu00=0.725, sigma2_2=0.148, phi2=23.5, tau2=0.349*0.148),
                   preferentiality.par = c(beta=-0.114))

# Creation of the mesh
mesh <- inla.mesh.2d(loc=galicia[,c("x","y")],
                     max.edge=c(25,30),cutoff = 5,
                     boundary=inla.mesh.segment(galicia.boundary))
par(mfrow=c(1,1))
plot(mesh,asp=1)

# Set of the number of simulations and burnin
control.mcmc <- control.mcmc.MCML(n.sim=12000, burnin=10000)

# A 3 by 3 km regular grid over Galicia, used to approximate the
# intractable integral of the log-Gaussian Cox process model
library(splancs)
grid.galicia <- gridpts(as.matrix(galicia.boundary),xs=3,ys=3)
points(grid.galicia,pch=20,cex=0.5)


# Fitting of the model
fit.ps <- lm.ps.MCML(log(lead)~1|1,coords=~x+y,data.response = galicia,
                     par0=par0,control.mcmc=control.mcmc,
                     which.is.preferential = (galicia$survey==1997)*1,
                     kappa1=1,kappa2=0.5,mesh=mesh,grid.intensity=grid.galicia)


# Summary of the model
summary(fit.ps)

# Spatial Prediction (on the log-scale)
pred.ps <- spatial.pred.lm.ps(fit.ps,target=3,control.mcmc=control.mcmc,standard.errors = TRUE,
                              quantiles=c(0.025,0.975),return.samples=TRUE)

# Predictive samples lead concentration
lead1997.samples <- exp(pred.ps$response$samples$preferential)
lead2000.samples <- exp(pred.ps$response$samples$non.preferential)

# Plotting of the predictions and standard errors
pred.ps$response$predictions$preferential <- apply(lead1997.samples,2,mean)
pred.ps$response$predictions$non.preferential <- apply(lead2000.samples,2,mean)
pred.ps$response$standard.errors$preferential <- apply(lead1997.samples,2,sd)
pred.ps$response$standard.errors$non.preferential <- apply(lead2000.samples,2,sd)

# To call the Help page on the plotting function below enter ?plot.pred.PrevMap.ps

# Predicted lead concentration in 1997 and 2000
plot(pred.ps,target=1,summary="predictions",main=c(1997,2000),zlim=c(1,8))

# Stndard errors for lead concentration in 1997 and 2000
plot(pred.ps,target=1,summary="standard.errors",main=c(1997,2000),zlim=c(0,2.55))

plot(pred.ps,target=2,summary="predictions",main=c("Sampling intensity 1997 survey (predictions)"))
plot(pred.ps,target=2,summary="standard.errors",main=c("Sampling intensity 1997 survey (standard errors)"))


