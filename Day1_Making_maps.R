rm(list=ls())
library(sf)
library(tmap)

rb <- read.csv("LiberiaRemoData.csv")

rb$prev <- rb$npos/rb$ntest

rb.sf <- st_as_sf(rb,coords = c("utm_x","utm_y"))
st_crs(rb.sf) <- 32629

liberia.adm0 <- st_read("./LBR_adm/LBR_adm0.shp")
liberia.adm0 <- st_transform(liberia.adm0,crs=32629)

map0 <- tm_shape(liberia.adm0) + 
  tm_borders(lwd=3) 
map0

map0+tm_shape(rb.sf)+tm_dots(size=0.5)

Map.with.points <- map0+tm_shape(rb.sf)+
  tm_bubbles("prev", col = "prev", 
             border.col = "black",
             style="fixed", 
             breaks=seq(0,0.4,0.05),
             palette="-RdYlBu",
             title.size="Prevalence", 
             scale = 1,
             title.col="") 
Map.with.points+
  tm_compass(type="8star", 
             position = c("right","top"))+
  tm_scale_bar(breaks = c(0,100,200),text.size=1,
               position=c("center","bottom")) 

tmap_mode("view")

Map.with.points
tmap_mode("plot")

liberia.wl <- st_read("./LBR_wat/LBR_water_lines_dcw.shp")
liberia.wl <- st_transform(liberia.wl,crs = 32629)

Map.with.points+tm_shape(liberia.wl)+
  tm_lines(col="blue", palette="dodgerblue3", 
           title.col="Waterways")

####
library(raster)

liberia.alt <- raster("./LBR_alt/LBR_alt.gri")
liberia.alt <- projectRaster(liberia.alt,
               crs=CRS("+init=epsg:32629"))

tm_shape(liberia.alt)+
tm_raster(title="Elevation")+
tm_shape(liberia.adm0)+
tm_borders(lwd=2)

liberia.alt <- mask(liberia.alt,
                    as(liberia.adm0,"Spatial"))

tm_shape(liberia.alt)+
tm_raster(title="Elevation (m)")+
tm_shape(liberia.adm0)+
tm_borders(lwd=2)

####
liberia.grid <- st_make_grid(liberia.adm0,
                             cellsize = 2000,
                             what="centers")

liberia.inout <- st_intersects(liberia.grid,
                               liberia.adm0,
                               sparse = FALSE)
liberia.grid <- liberia.grid[liberia.inout]

dist <- apply(st_distance(liberia.grid,
                          liberia.wl),1,min)/1000

dist.raster <- rasterFromXYZ(cbind(
  st_coordinates(liberia.grid),
  dist),
  crs="+init=epsg:32629")

tm_shape(dist.raster)+
  tm_raster(title="Distance from \n closest waterway (km)")+
  tm_shape(liberia.adm0)+
  tm_borders(lwd=2)+
  tm_shape(liberia.wl)+
  tm_lines(col="blue",lwd=2)

writeRaster(dist.raster,
            "Liberia_dist_cw.tiff",
            "GTiff")

liberia.max2km.dist <- 
extract(dist.raster,
        st_coordinates(liberia.grid),
        buffer=5000,
        fun=max)

dist.raster5 <- rasterFromXYZ(cbind(
  st_coordinates(liberia.grid),
  liberia.max5km.dist),
  crs="+init=epsg:32629")
plot(dist.raster5)

###
liberia.adm2 <- st_read("./LBR_adm/LBR_adm2.shp")
liberia.adm2 <- st_transform(liberia.adm2,crs=32629)

names.adm2 <- liberia.adm2$NAME_2
n.adm2 <- length(names.adm2)
mean.elev <- rep(NA,n.adm2)

liberia.adm2$Mean_elevation <- NA

for(i in 1:n.adm2) {
  ind.sel <- names.adm2==names.adm2[i]
  elev.r.i <- mask(liberia.alt,
                    as(
                    liberia.adm2[ind.sel,],
                    "Spatial"))
  liberia.adm2$Mean_elevation[ind.sel] <- 
    mean(values(elev.r.i),na.rm=TRUE)
}


map2 <- tm_shape(liberia.adm2) + 
  tm_borders(lwd=1) 

map2+tm_fill("Mean_elevation",
             title = "Mean elevation (m)")
