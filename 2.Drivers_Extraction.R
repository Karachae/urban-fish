# R version 4.1.2
# Extracting and re-projecting urbanization metrics
# Re-project to Lambert Cylindrical Equal Area

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ncdf4)
library(ncdf4.helpers)
library(tidyverse)
library(raster)
library(sp)
library(sf)
library(hierfstat)
library(rgdal)

d <- read.csv("LoopedGlobalMasterPopSheet.csv")
coordinates(d) <- ~ long + lat
proj4string(d)<- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84")
d<- spTransform(d, CRS("+proj=cea +lat_ts=45 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
st_crs(d)

## Human Population Density (2010) ##
## https://sedac.ciesin.columbia.edu/

# Human Population Density (2010) 25km buffer
raw_hupopdens_data <- "gpw_v4_population_density_rev11_2010_30_sec.tif"
hupopdens_data  <- raster(raw_hupopdens_data)
hupopdens_data 
projection(hupopdens_data) #check CRS
hupopdens_data <- projectRaster(hupopdens_data, crs="+proj=cea +lat_ts=45 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
plot(hupopdens_data)
ext1 <- extract(hupopdens_data,d,fun=mean,buffer=25000, na.rm=TRUE, df=TRUE) 
ext1
d_hupopdens25 <- cbind(d, ext1[,2])
names(d_hupopdens25)[ncol(d_hupopdens25)] <- "hupopdens25km"
d_hupopdens25 <- as.data.frame(d_hupopdens25)

# Human Population Density (2010) 50km buffer
ext2 <- extract(hupopdens_data,d,fun=mean,buffer=50000, na.rm=TRUE, df=TRUE) 
ext2
d_hupopdens50 <- cbind(d_hupopdens25, ext2[,2])
names(d_hupopdens50)[ncol(d_hupopdens50)] <- "hupopdens50km"
d_hupopdens50 <- as.data.frame(d_hupopdens50)

# Human Population Density (2010) 100km buffer
ext3 <- extract(hupopdens_data,d,fun=mean,buffer=100000, na.rm=TRUE, df=TRUE) 
ext3
d_hupopdens100 <- cbind(d_hupopdens50, ext3[,2])
names(d_hupopdens100)[ncol(d_hupopdens100)] <- "hupopdens100km"
d_hupopdens100 <- as.data.frame(d_hupopdens100)

# Human Population Density (2010) 200km buffer
ext4 <- extract(hupopdens_data,d,fun=mean,buffer=200000, na.rm=TRUE, df=TRUE) 
ext4
d_hupopdens200 <- cbind(d_hupopdens100, ext4[,2])
names(d_hupopdens200)[ncol(d_hupopdens200)] <- "hupopdens200km"
d_hupopdens200 <- as.data.frame(d_hupopdens200)

## Halpern et al. 2015 ##
## https://doi.org/10.1038/ncomms8615 

# Re-projected on R server (version 3.6.3) because of memory issues
raw_CumImp_data <- "global_cumul_impact_2013_all_layers.tif"
CumImp_data <- raster(raw_CumImp_data)
CumImp_data
projection(CumImp_data) #check CRS
CumImp_data <- projectRaster(CumImp_data, crs="+proj=cea +lat_ts=45 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
writeRaster(CumImp_data, "global_cumul_impact_2013_all_layers_reproj", format="GTiff")

# Cumulative Impacts 25km buffer
raw_CumImp_data <- "global_cumul_impact_2013_all_layers_reproj.tif"
CumImp_data <- raster(raw_CumImp_data)
CumImp_data
projection(CumImp_data) #check CRS
plot(CumImp_data)
ext5 <- extract(CumImp_data,d,fun=mean,buffer=25000, na.rm=TRUE, df=TRUE) 
ext5
d_hupopdens.CI25 <- cbind(d_hupopdens200, ext5[,2])
names(d_hupopdens.CI25)[ncol(d_hupopdens.CI25)] <- "impacts25km"
d_hupopdens.CI25 <- as.data.frame(d_hupopdens.CI25)

# Cumulative Impacts 50km buffer
ext6 <- extract(CumImp_data,d,fun=mean,buffer=50000, na.rm=TRUE, df=TRUE) 
ext6
d_hupopdens.CI50 <- cbind(d_hupopdens.CI25, ext6[,2])
names(d_hupopdens.CI50)[ncol(d_hupopdens.CI50)] <- "impacts50km"
d_hupopdens.CI50 <- as.data.frame(d_hupopdens.CI50)

# Cumulative Impacts 100km buffer
ext7 <- extract(CumImp_data,d,fun=mean,buffer=100000, na.rm=TRUE, df=TRUE) 
ext7
d_hupopdens.CI100 <- cbind(d_hupopdens.CI50, ext7[,2])
names(d_hupopdens.CI100)[ncol(d_hupopdens.CI100)] <- "impacts100km"
d_hupopdens.CI100 <- as.data.frame(d_hupopdens.CI100)

# Cumulative Impacts 200km buffer
ext8 <- extract(CumImp_data,d,fun=mean,buffer=200000, na.rm=TRUE, df=TRUE) 
ext8
d_hupopdens.CI200 <- cbind(d_hupopdens.CI100, ext8[,2])
names(d_hupopdens.CI200)[ncol(d_hupopdens.CI200)] <- "impacts200km"
d_hupopdens.CI200 <- as.data.frame(d_hupopdens.CI200)

## Kroodsma et al. 2016 ##
## DOI: 10.1126/science.aao5646

# Fishing pressure (hours of fishing per km2) 25km buffer
raw_fishing_data <- "fishing_raster_05_20180228.tif"
fishing_data<- raster(raw_fishing_data)
fishing_data
projection(fishing_data) #check CRS
st_crs(fishing_data) # CRS already in format we want
plot(fishing_data)
ext9 <- extract(fishing_data,d,fun=mean,buffer=25000, na.rm=TRUE, df=TRUE) 
ext9
d_hupopdens.CI.fishing25 <- cbind(d_hupopdens.CI200, ext9[,2])
names(d_hupopdens.CI.fishing25)[ncol(d_hupopdens.CI.fishing25)] <- "fishing_effort25km"
d_hupopdens.CI.fishing25 <- as.data.frame(d_hupopdens.CI.fishing25)

# Fishing pressure (hours of fishing per km2) 50km buffer
ext10 <- extract(fishing_data,d,fun=mean,buffer=50000, na.rm=TRUE, df=TRUE) 
ext10
d_hupopdens.CI.fishing50 <- cbind(d_hupopdens.CI.fishing25, ext10[,2])
names(d_hupopdens.CI.fishing50 )[ncol(d_hupopdens.CI.fishing50 )] <- "fishing_effort50km"
d_hupopdens.CI.fishing50  <- as.data.frame(d_hupopdens.CI.fishing50)

# Fishing pressure (hours of fishing per km2) 100km buffer
ext11 <- extract(fishing_data,d,fun=mean,buffer=100000, na.rm=TRUE, df=TRUE) 
ext11
d_hupopdens.CI.fishing100 <- cbind(d_hupopdens.CI.fishing50, ext11[,2])
names(d_hupopdens.CI.fishing100)[ncol(d_hupopdens.CI.fishing100)] <- "fishing_effort100km"
d_hupopdens.CI.fishing100 <- as.data.frame(d_hupopdens.CI.fishing100)

# Fishing pressure (hours of fishing per km2) 200km buffer
ext12 <- extract(fishing_data,d,fun=mean,buffer=200000, na.rm=TRUE, df=TRUE) 
ext12
d_hupopdens.CI.fishing200 <- cbind(d_hupopdens.CI.fishing100, ext12[,2])
names(d_hupopdens.CI.fishing200)[ncol(d_hupopdens.CI.fishing200)] <- "fishing_effort200km"
d_hupopdens.CI.fishing200 <- as.data.frame(d_hupopdens.CI.fishing200)

# Check everything is in same CRS
st_crs(d) == st_crs(CumImp_data) 
st_crs(CumImp_data) == st_crs(hupopdens_data) 
st_crs(hupopdens_data) == st_crs(fishing_data) 
st_crs(fishing_data) == st_crs(CumImp_data) 

### Final Dataset ### 

write.csv(d_hupopdens.CI.fishing200, "ReProj_FinalUrbanizationdata.csv", row.names = FALSE)
