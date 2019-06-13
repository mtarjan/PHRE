##PHRE in hyenas
##Prepared by M Tarjan
##June 11, 2019

##INPUTS
##locs: dataframe with colnames x & y
##rast: raster with 1 or more layers
##smoother: bandwidth value (matric or scalar) or 'default'
##percent: percent kernel to calculate (e.g., 90)

##prepare inputs
data<-read.csv("C:/Users/max/Desktop/Tarjan/hyena_data/OHB09.csv")
plot(data$LON, data$LAT)

fence<-rgdal::readOGR(dsn = "C:/Users/max/Desktop/Tarjan/hyena_data/Fences/Fences", layer="Fences")

##reproject to projected coordinate system
##original - EPSG:4326 WGS 84.  convert- EPSG:32733 UTM 33S
raster::crs(fence)
fence.utm<-spTransform(fence, crs(rst))
