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
plot(data$LAT, data$LON)

fence<-rgdal::readOGR(dsn = "C:/Users/max/Desktop/Tarjan/hyena_data", layer="ONR_shape")
