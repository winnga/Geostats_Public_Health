rm(list=ls())
library(sf)
library(tmap)
rb <- read.csv("LiberiaRemoData.csv")

rb$Prevalence <- rb$npos/rb$ntest

Liberia.adm0 <- st_read("LBR_adm/LBR_adm0.shp")

Map.adm0.border <- 
  tm_shape(Liberia.adm0) + 
  tm_borders() 


rb.sf <- st_as_sf(rb, coords = c("long","lat"),
                  crs=st_crs(Liberia.adm0)$proj4string)

# Simple points map
Map.adm0.border+tm_shape(rb.sf)+
  tm_dots(size=1)

Map.with.points <- Map.adm0.border+tm_shape(rb.sf)+
  tm_bubbles("Prevalence", col = "Prevalence", 
           border.col = "black",
           style="fixed", 
           breaks=seq(0,0.4,0.05),
           palette="-RdYlBu",
           title.size="Prevalence", 
           scale = 1.5,
           title.col="") 

Map.with.points+  
  tm_compass(type="8star", 
             position = c("right","top"))+
  tm_scale_bar(breaks = c(0,100,200),size=1,
               position=c("center","bottom")) 

tmap_mode("view")
Map.with.points

Liberia.grid.sq <- st_make_grid(Liberia.adm0,
                                cellsize = 0.07,
                                what="centers")

Liberia.inout <- st_intersects(Liberia.grid.sq,
                               Liberia.adm0,
                               sparse = FALSE)
Liberia.grid <- Liberia.grid.sq[Liberia.inout]

water <- st_read("LBR_wat/LBR_water_lines_dcw.shp")

dist <- apply(st_distance(Liberia.grid,water),1,min)/1000
library(raster)
dist.raster <- rasterFromXYZ(cbind(
  st_coordinates(Liberia.grid),
  dist),
  crs="+init=epsg:4326")

tm_shape(dist.raster)+tm_raster(title="Distance from river (km)")+
  tm_shape(water)+tm_lines(col="blue")

rb$dist_water <- extract(dist.raster,rb.sf)

plot(rb$dist_water,rb$Prevalence,
     col=1*(rb$elevation>50)+2)


### Repeat using a UTM projection

Liberia.adm0.utm <- st_transform(Liberia.adm0,32629)
st_crs(Liberia.adm0.utm)
