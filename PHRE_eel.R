##apply phre function to eel data
##M Tarjan
##July 11, 2019

library(rgeos)
library(readxl)
library(ks)
##load data; specify file paths to supporting information 9 data
data <- read.csv2("C:/Users/max/Downloads/hr_33.csv")
data<-cbind(data$x, data$y)
grid<-read.csv2("C:/Users/max/Downloads/habitat.csv")

##create list of rasters from the grid array
library(raster)
rasters<-list(0)
rasters[[1]]<-rasterFromXYZ(cbind(grid$x, grid$y, grid$RivLength))
rasters[[2]]<-rasterFromXYZ(cbind(grid$x, grid$y, log(grid$BankDist+0.01)))

##specify the baseline smoothing parameters
#smoother<-cbind(c(1,0), c(0,1)) #(ATOS, 0) (0, log(dist))

##apply phre function
HR<-phre(locs=data, rast=rasters, smooth='default', percent=90, resolution = 1000)

##plot of phre list objects
##zoomed in on polygons
plot(HR$Poly); plot(HR$array, add=T); plot(HR$Polygon, col="transparent", border='red', add=T); points(HR$locs, pch=20, cex=.1, col="black")
