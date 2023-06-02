rm(list=ls())
library(PrevMap)
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
