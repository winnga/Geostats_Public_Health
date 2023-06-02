rm(list=ls())

rb <- read.csv("LiberiaRemoData.csv")


glm.fit1 <- glm(cbind(npos,ntest-npos) ~ I(utm_x/1000)+
                                         I(utm_y/1000),
                data=rb, family=binomial)

summary(glm.fit1)

library(sf)
liberia.adm0 <- st_read("./LBR_adm/LBR_adm0.shp")
liberia.adm0 <- st_transform(liberia.adm0,32629)

library(splancs)
liberia.grid <- 
gridpts(st_coordinates(liberia.adm0)[,1:2],
        xs=3000,ys=3000)

predict.prev1 <- 
predict(glm.fit1,newdata=data.frame(utm_x=liberia.grid[,1],
                                    utm_y=liberia.grid[,2]),
        type="response")
r.prev1 <- rasterFromXYZ(cbind(liberia.grid,Prevalence=predict.prev1),
                         crs=CRS("+init=epsg:32629"))

library(tmap)
tm_shape(r.prev1)+tm_raster(col="Prevalence",
                            breaks=seq(0,0.35,length=10))

beta.hat <- coef(glm.fit1)
Sigma.par <- vcov(glm.fit1)
Sigma.par.sroot <- t(chol(Sigma.par))
n.sim <- 1000
beta.samples <- sapply(1:n.sim,
                       function(i)
                         beta.hat+Sigma.par.sroot%*%rnorm(3))

prev.samples <- sapply(1:n.sim, function(i) {
                                  lp.sample.i <- beta.samples[1,i]+
                                  beta.samples[2,i]*liberia.grid[,1]/1000+
                                  beta.samples[3,i]*liberia.grid[,2]/1000
                                  
                                  out <- exp(lp.sample.i)/(1+exp(lp.sample.i))
                                  return(out)
                                })

prob.exceed20.1 <- apply(prev.samples,1,function(x) mean(x>0.2))
r.exceed20.1 <- rasterFromXYZ(cbind(liberia.grid,
                                    Exceed.20=prob.exceed20.1),
                         crs=CRS("+init=epsg:32629"))
tm_shape(r.exceed20.1)+tm_raster(col="Exceed.20")

###
glm.fit2 <- glm(cbind(npos,ntest-npos)~log(elevation),
               data=rb,family=binomial)

elevation <- raster("./LBR_alt/LBR_alt.gri")
elevation <- projectRaster(elevation,crs=CRS("+init=epsg:32629"))

elevation.pred <- extract(elevation,liberia.grid)
ind.na <- which(is.na(elevation.pred))
elevation.pred <- elevation.pred[-ind.na]
liberia.grid2 <- liberia.grid[-ind.na,]

predict.prev2 <- 
  predict(glm.fit2,newdata=data.frame(elevation=elevation.pred),
          type="response")
r.prev2 <- rasterFromXYZ(cbind(liberia.grid2,Prevalence=predict.prev2),
                         crs=CRS("+init=epsg:32629"))

library(tmap)
tm_shape(r.prev2)+tm_raster(col="Prevalence",
                            breaks=seq(0,0.35,length=10))

beta.hat <- coef(glm.fit2)
Sigma.par <- vcov(glm.fit2)
Sigma.par.sroot <- t(chol(Sigma.par))
n.sim <- 1000
beta.samples <- sapply(1:n.sim,
                       function(i)
                         beta.hat+Sigma.par.sroot%*%rnorm(2))
prev.samples <- sapply(1:n.sim, function(i) {
  lp.sample.i <- beta.samples[1,i]+
    beta.samples[2,i]*log(elevation.pred)
  
  out <- exp(lp.sample.i)/(1+exp(lp.sample.i))
  return(out)
})


prob.exceed20.2 <- apply(prev.samples,1,function(x) mean(x>0.2))
r.exceed20.2 <- rasterFromXYZ(cbind(liberia.grid2,
                                    Exceed.20=prob.exceed20.2),
                              crs=CRS("+init=epsg:32629"))
tm_shape(r.exceed20.2)+tm_raster(col="Exceed.20")

rb$prev <- rb$npos/rb$ntest
rb$prev.hat1 <- extract(r.prev1,rb[,c("utm_x","utm_y")])
rb$prev.hat2 <- extract(r.prev2,rb[,c("utm_x","utm_y")])

plot(rb$prev.hat1,rb$prev)
abline(0,1)

plot(rb$prev.hat2,rb$prev)
abline(0,1)

glm.fit1$aic
glm.fit2$aic

liberia.adm2 <- st_read("./LBR_adm/LBR_adm2.shp")
liberia.adm2 <- st_transform(liberia.adm2,32629)

liberia.grid.sf <- st_as_sf(data.frame(
                   utm_x=liberia.grid2[,1],
                   utm_y=liberia.grid2[,2]),
                   coords=c("utm_x","utm_y"))
st_crs(liberia.grid.sf) <- 32629

prev.hat.distr <- rep(NA,nrow(liberia.adm2))
prev.q025.distr <- rep(NA,nrow(liberia.adm2))
prev.q975.distr <- rep(NA,nrow(liberia.adm2))

for(i in 1:nrow(liberia.adm2)) {
  inter.i <- st_intersects(liberia.adm2[i,],liberia.grid.sf)
  prev.distr.samples.i <- apply(prev.samples[inter.i[[1]],],2,mean)
  prev.hat.distr[i] <- mean(prev.distr.samples.i)
  prev.q025.distr[i] <- quantile(prev.distr.samples.i,0.025)
  prev.q025.distr[i] <- quantile(prev.distr.samples.i,0.975)
}
 
liberia.adm2$Prevalence <- prev.hat.distr

tm_shape(liberia.adm2)+tm_borders()+
  tm_fill(col="Prevalence")
###

rb <- read.csv("LiberiaRemoData.csv")
rb$ID <- 1:nrow(rb)
glmer.fit <- glmer(cbind(npos,ntest-npos) ~ log(elevation)+(1|ID),
                   data=rb,family=binomial)

Z.hat <- ranef(glmer.fit)$ID[,1]

coords <- rb[,c("utm_x","utm_y")]

data.variogram <- data.frame(Z.hat=Z.hat,
                             utm_x=coords[,1],
                             utm_y=coords[,2])
vari <- variogram(data.variogram,~Z.hat,
                  ~I(utm_x/1000)+I(utm_y/1000))
plot(vari,xlab="Distance (km)",ylab="Variogram",
     type="b")

spat.corr.diagnostic(npos~log(elevation),units.m=~ntest,
                     coords=~I(utm_x/1000)+I(utm_y/1000),
                     data=rb,likelihood = "Binomial",
                     which.test = "variogram")
