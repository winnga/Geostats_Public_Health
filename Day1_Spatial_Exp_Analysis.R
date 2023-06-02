rm(list=ls())
library(PrevMap)
maln <- read.csv("malnutrition.csv")

maln$ID.clusters <- create.ID.coords(maln,~utm_x+utm_y)

library(lme4)

plot(HAZ ~ age, data=maln,
     xlab="Age (years)", ylab="HAZ")
max.vec <- function(vec,c) 
           sapply(vec,function(x) max(x,c))

lmer.fit <- 
lmer(HAZ ~ age + max.vec(age,2) + (1|ID.clusters), 
     data = maln)

summary(lmer.fit)

Z.hat <- ranef(lmer.fit)$ID.clusters[,1]
coords <- unique(maln[,c("utm_x","utm_y")])
data.variogram <- data.frame(Z.hat=Z.hat,
                             utm_x=coords[,1],
                             utm_y=coords[,2])

variogram.Z <- 
variogram(data.variogram,
          ~Z.hat,~utm_x+utm_y,
          uvec=seq(0,300,length=15))
plot(variogram.Z,typ="b")

spat.corr.diagnostic(HAZ ~ age + max.vec(age,2),
                     coords=~utm_x+utm_y,
                     data=maln,
                     ID.coords = maln$ID.clusters,
                     likelihood = "Gaussian")
###
mosq <- read.csv("anopheles.csv")

par(mfrow=c(1,1),mar=c(4,4,2,1))
elev <- read.csv("anopheles_elevation.csv")
elev <- rasterFromXYZ(elev)
plot(elev)
points(mosq[,c("web_x","web_y")],pch=20)
plot(log(An.gambiae) ~ elevation, data = mosq,pch=20,cex=0.5)
glm.fit <- glm(An.gambiae ~ elevation, data = mosq, family = poisson)
summary(glm.fit)
abline(glm.fit)

check.spat <- spat.corr.diagnostic(An.gambiae~elevation,
                                   coords = ~I(web_x/1000)+I(web_y/1000),
                                   data=mosq,likelihood = "Poisson",
                                   uvec=seq(0,10,length=15),n.sim=10000)
