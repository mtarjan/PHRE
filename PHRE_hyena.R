##PHRE in hyenas
##Prepared by M Tarjan
##June 11, 2019

##INPUTS
##locs: dataframe with colnames x & y
##rast: raster with 1 or more layers
##smoother: bandwidth value (matric or scalar) or 'default'
##percent: percent kernel to calculate (e.g., 90)

library(maptools) ##required for spatialpoints
library(raster)

##prepare inputs
data<-read.csv("C:/Users/max/Desktop/Tarjan/hyena_data/OHB08.csv")
#plot(data$LON, data$LAT)

locs.ll<-subset(data, select=c("LON", "LAT"))
locs.ll<-SpatialPoints(locs.ll)
crs(locs.ll)<-CRS("+init=epsg:4326") ##assign a coordinate system
locs.utm<-spTransform(locs.ll, CRSobj = CRS("+init=epsg:32733")) ##reproject

fence<-rgdal::readOGR(dsn = "C:/Users/max/Desktop/Tarjan/hyena_data/Fences/Fences", layer="Fences")
out.fence.utm<-rgdal::readOGR(dsn = "C:/Users/max/Desktop/Tarjan/hyena_data/Fences/Fences", layer="outside_fences_utm2")

##reproject to projected coordinate system
##original - EPSG:4326 WGS 84.  convert- EPSG:32733 UTM 33S
raster::crs(fence)
fence.utm<-spTransform(fence, CRSobj = CRS("+init=epsg:32733"))
out.fence.utm<-spTransform(out.fence.utm, CRSobj = CRS("+init=epsg:32733"))

plot(out.fence.utm); plot(locs.utm, add=T)

##create habitat rasters
bbox<-summary(out.fence.utm)$bbox
habitat<-raster(xmn = bbox[1,1]+500, xmx = bbox[1,2]-500, ymn = bbox[2,1]+500, ymx = bbox[2,2]-500, crs = crs(out.fence.utm), resolution = 50)
dd<-rgeos::gDistance(spgeom1=out.fence.utm, spgeom2=as(habitat, "SpatialPoints"), byid=T) ##calc distance from fence for every point in raster
habitat[] = log(apply(dd,1,min)+0.01)
plot(habitat); plot(fence.utm, add=T); plot(locs.utm, add=T)

#plot(HR$Poly); plot(habitat, add=T); plot(fence.utm, add=T)

rasters<-list(x=habitat, y=habitat, z=habitat)
rasters[[1]][]<-rasterToPoints(habitat)[,1]
rasters[[2]][]<-rasterToPoints(habitat)[,2]
#rasters[[3]]<-habitat

rm(habitat, dd) ##remove habitat from workspace to save space

##run PHRE
##apply phre function
HR<-phre(locs=locs.utm@coords, rast=rasters, smooth='default', percent=90, resolution=1500)
#HR<-phre(locs=locs.utm@coords, rast=rasters)##resolution here is equal to the number of cells across the input raster, so a higher number leads to better resolution; default is 500

#gc() ##garbage collection to free up memory

##for trouble-shooting in phre function
#locs<-locs.utm@coords; rast<-rasters; smooth<-'default'; percent=90; resolution=250

##plot of phre list objects
##zoomed in on home range
plot(HR$Poly); plot(HR$array, add=T); plot(HR$Polygon, col="transparent", border='red', add=T); points(HR$locs, pch=20, cex=.1, col="black"); plot(fence.utm, add=T)

##zoom on park
plot(fence.utm); plot(HR$array, add=T); plot(HR$Polygon, col="transparent", border='red', add=T); points(HR$locs, pch=20, cex=.1, col="black"); plot(fence.utm, add=T)
