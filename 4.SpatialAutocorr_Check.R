# R version 3.6.3 (on server)
# Testing for spatial autocorrelation

library(ggplot2)
library(marmap)
library(dplyr)
library(spdep)

## Calculate Distance Matrix ##

# download ocean depth map
ocean_map <- getNOAA.bathy(lon1 = -180, lon2 = 180, lat1 = 90, lat2 = -90, resolution = 4)
plot(ocean_map)

#format map
ocean_map <- marmap::as.raster(ocean_map)
ocean_map[ocean_map < 0] <- NA
ocean_map[is.na(ocean_map)] <- -500
ocean_map[ocean_map > 0] <- NA
ocean_map[is.na(ocean_map)] <- 500
ocean_map <- as.bathy(ocean_map)
plot(ocean_map)

# upload/import data/load site info
dataset<-  read.csv("LoopedGlobalMasterPopSheet.csv")
sites <- dplyr::select(dataset, pop, lat, long)

# plot with ggplot2
nbw_map_ocean <- autoplot(ocean_map, geom=c("raster")) +
  geom_point(data=sites, aes(x=long, y=lat), pch=8, size=1, stroke=0.5, colour="yellow")
nbw_map_ocean

# prepare site coords in necessary format
site_coords <- sites[,c("long", "lat")]
colnames(site_coords) <- c("x", "y")

# calculate distance
trans <- trans.mat(ocean_map, min.depth = 0)
# actual distance
dist <- lc.dist(trans,site_coords,res="dist") # run overnight 
max(dist) # check for unrealistic values

## Deal With NAs in FST & Ne

xy <- dplyr::select(dataset, long, lat) # coordinates
DistTable <- as.matrix(dist)

#find lines that have NA in FST
Fst <- dataset$global_fst
Fst
om_fst <- which(is.na(Fst), arr.ind=TRUE)
#remove from the distance matrix
dist
DistTable1 <- DistTable[-(om_fst), -(om_fst)]
dist1 <- as.dist(DistTable1)

#find lines that have NA in Ne
Ne <-dataset$Ne
Ne
om_ne <-which(is.na(Ne), arr.ind=TRUE)
#remove from the distance matrix
DistTable2 <- DistTable[-(om_ne), -(om_ne)]
dist2 <- as.dist(DistTable2)

## Convert to Appropriate Format ##

# Gene Diversity & Allelic Richness (because no NAs in either)
dist_m<- as.matrix(dist) 
dist_lw <-mat2listw(dist_m)
# FST
dist_m1<- as.matrix(dist1) 
dist_lw1 <-mat2listw(dist_m1)
# Ne
dist_m2<- as.matrix(dist2) 
dist_lw2 <-mat2listw(dist_m2)

## Check for spatial autocorrelation ##

# Gene Diversity Models No Priors #

# Human Population Density
res.MH_humanpopdens200 <- residuals(MH_humanpopdens200)[,1]
# removing rows from distance matrix because NAs in data
om_hupopdens200 <-which(is.na(fishy.he$hupopdens200km), arr.ind=TRUE)
dist_m_hupopdens200 <- dist_m[-(om_hupopdens200),-(om_hupopdens200)]
dist_lw_hupopdens200 <- mat2listw(dist_m_hupopdens200)
moran.test(res.MH_humanpopdens200, dist_lw_hupopdens200) # no, I = -3.104471e-03

res.MH_humanpopdens100 <- residuals(MH_humanpopdens100)[,1]
# removing rows from distance matrix because NAs in data
om_hupopdens100 <-which(is.na(fishy.he$hupopdens100km), arr.ind=TRUE)
dist_m_hupopdens100 <- dist_m[-(om_hupopdens100),-(om_hupopdens100)]
dist_lw_hupopdens100 <- mat2listw(dist_m_hupopdens100)
moran.test(res.MH_humanpopdens100, dist_lw_hupopdens100) # no, I = -3.618144e-03 

res.MH_humanpopdens50 <- residuals(MH_humanpopdens50)[,1]
# removing rows from distance matrix because NAs in data
om_hupopdens50 <-which(is.na(fishy.he$hupopdens50km), arr.ind=TRUE)
dist_m_hupopdens50 <- dist_m[-(om_hupopdens50),-(om_hupopdens50)]
dist_lw_hupopdens50 <- mat2listw(dist_m_hupopdens50)
moran.test(res.MH_humanpopdens50, dist_lw_hupopdens50) # no, I = -4.839801e-03

res.MH_humanpopdens25 <- residuals(MH_humanpopdens25)[,1]
# removing rows from distance matrix because NAs in data
om_hupopdens25 <-which(is.na(fishy.he$hupopdens25km), arr.ind=TRUE)
dist_m_hupopdens25 <- dist_m[-(om_hupopdens25),-(om_hupopdens25)]
dist_lw_hupopdens25 <- mat2listw(dist_m_hupopdens25)
moran.test(res.MH_humanpopdens25, dist_lw_hupopdens25) # no, I = -4.970371e-03

# Cumulative Human Impacts
res.MH_huimpacts200 <- residuals(MH_huimpacts200)[,1]
moran.test(res.MH_huimpacts200, dist_lw) # no, I = -4.129590e-03

res.MH_huimpacts100 <- residuals(MH_huimpacts100)[,1]
moran.test(res.MH_huimpacts100, dist_lw) # no, I = -4.388479e-03

res.MH_huimpacts50 <- residuals(MH_huimpacts50)[,1]
moran.test(res.MH_huimpacts50, dist_lw) # no, I = -4.699615e-03

res.MH_huimpacts25 <- residuals(MH_huimpacts25)[,1]
moran.test(res.MH_huimpacts25, dist_lw) # no, I = -4.750789e-03

# Fishing Pressure
res.MH_fishing200 <- residuals(MH_fishing200)[,1]
# removing rows from distance matrix because NAs in data
om_fishing200 <-which(is.na(fishy.he$fishing_effort200km), arr.ind=TRUE)
dist_m_fishing200 <- dist_m[-(om_fishing200),-(om_fishing200)]
dist_lw_fishing200 <- mat2listw(dist_m_fishing200)
moran.test(res.MH_fishing200, dist_lw_fishing200) # no, I = -2.745815e-03

res.MH_fishing100 <- residuals(MH_fishing100)[,1]
# removing rows from distance matrix because NAs in data
om_fishing100 <-which(is.na(fishy.he$fishing_effort100km), arr.ind=TRUE)
dist_m_fishing100 <- dist_m[-(om_fishing100),-(om_fishing100)]
dist_lw_fishing100 <- mat2listw(dist_m_fishing100)
moran.test(res.MH_fishing100, dist_lw_fishing100) # no, I = -3.547316e-03

res.MH_fishing50 <- residuals(MH_fishing50)[,1]
# removing rows from distance matrix because NAs in data
om_fishing50 <-which(is.na(fishy.he$fishing_effort50km), arr.ind=TRUE)
dist_m_fishing50 <- dist_m[-(om_fishing50),-(om_fishing50)]
dist_lw_fishing50 <- mat2listw(dist_m_fishing50)
moran.test(res.MH_fishing50, dist_lw_fishing50) # no, I = -4.341077e-03

res.MH_fishing25 <- residuals(MH_fishing25)[,1]
# removing rows from distance matrix because NAs in data
om_fishing25 <-which(is.na(fishy.he$fishing_effort25km), arr.ind=TRUE)
dist_m_fishing25 <- dist_m[-(om_fishing25),-(om_fishing25)]
dist_lw_fishing25 <- mat2listw(dist_m_fishing25)
moran.test(res.MH_fishing25, dist_lw_fishing25) # no, I = -3.851053e-03

# Allelic Richness Models No Priors #

# Human Population Density
res.MAR_humanpopdens200 <- residuals(MAR_humanpopdens200)[,1]
moran.test(res.MAR_humanpopdens200, dist_lw_hupopdens200) # no, I = -3.147267e-03

res.MAR_humanpopdens100 <- residuals(MAR_humanpopdens100)[,1]
moran.test(res.MAR_humanpopdens100, dist_lw_hupopdens100) # no, I = -3.471933e-03

res.MAR_humanpopdens50 <- residuals(MAR_humanpopdens50)[,1]
moran.test(res.MAR_humanpopdens50, dist_lw_hupopdens50) # no, I = -4.819972e-03 

res.MAR_humanpopdens25 <- residuals(MAR_humanpopdens25)[,1]
moran.test(res.MAR_humanpopdens25, dist_lw_hupopdens25) # no, I =  -4.855876e-03 

# Cumulative Human Impacts
res.MAR_huimpacts200 <- residuals(MAR_huimpacts200)[,1]
moran.test(res.MAR_huimpacts200, dist_lw) # no, I =  -3.625336e-03

res.MAR_huimpacts100 <- residuals(MAR_huimpacts100)[,1]
moran.test(res.MAR_huimpacts100, dist_lw) # no, I = -4.229930e-03

res.MAR_huimpacts50 <- residuals(MAR_huimpacts50)[,1]
moran.test(res.MAR_huimpacts50, dist_lw) # no, I = -4.569538e-03

res.MAR_huimpacts25 <- residuals(MAR_huimpacts25)[,1]
moran.test(res.MAR_huimpacts25, dist_lw) # no, I = -4.697962e-03

# Fishing Pressure
res.MAR_fishing200 <- residuals(MAR_fishing200)[,1]
moran.test(res.MAR_fishing200, dist_lw_fishing200) # no, I = -2.803954e-03

res.MAR_fishing100 <- residuals(MAR_fishing100)[,1]
moran.test(res.MAR_fishing100, dist_lw_fishing100) # no, I = -3.786989e-03

res.MAR_fishing50 <- residuals(MAR_fishing50)[,1]
moran.test(res.MAR_fishing50, dist_lw_fishing50) # no, I = -4.736018e-03

res.MAR_fishing25 <- residuals(MAR_fishing25)[,1]
moran.test(res.MAR_fishing25, dist_lw_fishing25) # no, I = -2.622460e-03

# Fst Models No Priors #

# Human Population Density
res.MFST_humanpopdens200 <- residuals(MFST_humanpopdens200)[,1]
# removing rows from distance matrix because NAs in data
om1_hupopdens200 <-which(is.na(fishy.fst$hupopdens200km), arr.ind=TRUE)
dist_m1_hupopdens200 <- dist_m1[-(om1_hupopdens200), -(om1_hupopdens200)]
dist_lw1_hupopdens200 <- mat2listw(dist_m1_hupopdens200)
moran.test(res.MFST_humanpopdens200, dist_lw1_hupopdens200) # no, I = -7.222203e-04

res.MFST_humanpopdens100 <- residuals(MFST_humanpopdens100)[,1]
# removing rows from distance matrix because NAs in data
om1_hupopdens100 <-which(is.na(fishy.fst$hupopdens100km), arr.ind=TRUE)
dist_m1_hupopdens100 <- dist_m1[-(om1_hupopdens100), -(om1_hupopdens100)]
dist_lw1_hupopdens100 <- mat2listw(dist_m1_hupopdens100)
moran.test(res.MFST_humanpopdens100, dist_lw1_hupopdens100) # no, I = -9.140351e-04

res.MFST_humanpopdens50 <- residuals(MFST_humanpopdens50)[,1]
# removing rows from distance matrix because NAs in data
om1_hupopdens50 <-which(is.na(fishy.fst$hupopdens50km), arr.ind=TRUE)
dist_m1_hupopdens50 <- dist_m1[-(om1_hupopdens50), -(om1_hupopdens50)]
dist_lw1_hupopdens50 <- mat2listw(dist_m1_hupopdens50)
moran.test(res.MFST_humanpopdens50, dist_lw1_hupopdens50) # no, I = -9.980684e-04

res.MFST_humanpopdens25 <- residuals(MFST_humanpopdens25)[,1]
# removing rows from distance matrix because NAs in data
om1_hupopdens25 <-which(is.na(fishy.fst$hupopdens25km), arr.ind=TRUE)
dist_m1_hupopdens25 <- dist_m1[-(om1_hupopdens25), -(om1_hupopdens25)]
dist_lw1_hupopdens25 <- mat2listw(dist_m1_hupopdens25)
moran.test(res.MFST_humanpopdens25, dist_lw1_hupopdens25) # no, I = -1.065673e-03

# Cumulative Human Impacts
res.MFST_huimpacts200 <- residuals(MFST_huimpacts200)[,1]
moran.test(res.MFST_huimpacts200, dist_lw1) # no, I = -1.042395e-03

res.MFST_huimpacts100 <- residuals(MFST_huimpacts100)[,1]
moran.test(res.MFST_huimpacts100, dist_lw1) # no, I = -1.080554e-03

res.MFST_huimpacts50 <- residuals(MFST_huimpacts50)[,1]
moran.test(res.MFST_huimpacts50, dist_lw1) # no, I = -1.087225e-03

res.MFST_huimpacts25 <- residuals(MFST_huimpacts25)[,1]
moran.test(res.MFST_huimpacts25, dist_lw1) # no, I = -1.097912e-03

# Fishing Pressure
res.MFST_fishing200 <- residuals(MFST_fishing200)[,1]
# removing rows from distance matrix because NAs in data
om1_fishing200 <-which(is.na(fishy.fst$fishing_effort200km), arr.ind=TRUE)
dist_m1_fishing200 <- dist_m1[-(om1_fishing200),-(om1_fishing200)]
dist_lw1_fishing200 <- mat2listw(dist_m1_fishing200)
moran.test(res.MFST_fishing200, dist_lw1_fishing200) # no, I = -1.348302e-03

res.MFST_fishing100 <- residuals(MFST_fishing100)[,1]
# removing rows from distance matrix because NAs in data
om1_fishing100 <-which(is.na(fishy.fst$fishing_effort100km), arr.ind=TRUE)
dist_m1_fishing100 <- dist_m1[-(om1_fishing100),-(om1_fishing100)]
dist_lw1_fishing100 <- mat2listw(dist_m1_fishing100)
moran.test(res.MFST_fishing100, dist_lw1_fishing100) # no, I = -1.463013e-03 

res.MFST_fishing50 <- residuals(MFST_fishing50)[,1]
# removing rows from distance matrix because NAs in data
om1_fishing50 <-which(is.na(fishy.fst$fishing_effort50km), arr.ind=TRUE)
dist_m1_fishing50 <- dist_m1[-(om1_fishing50),-(om1_fishing50)]
dist_lw1_fishing50 <- mat2listw(dist_m1_fishing50)
moran.test(res.MFST_fishing50, dist_lw1_fishing50) # no, I = -1.262671e-03

res.MFST_fishing25 <- residuals(MFST_fishing25)[,1]
# removing rows from distance matrix because NAs in data
om1_fishing25 <-which(is.na(fishy.fst$fishing_effort25km), arr.ind=TRUE)
dist_m1_fishing25 <- dist_m1[-(om1_fishing25),-(om1_fishing25)]
dist_lw1_fishing25 <- mat2listw(dist_m1_fishing25)
moran.test(res.MFST_fishing25, dist_lw1_fishing25) # no, I = -8.844649e-04   

# Ne Models No Priors #
 
# Human Population Density
res.MNE_humanpopdens200 <- residuals(MNE_humanpopdens200)[,1]
# removing rows from distance matrix because NAs in data
om2_hupopdens200 <-which(is.na(fishy.ne$hupopdens200km), arr.ind=TRUE)
dist_m2_hupopdens200 <- dist_m2[-(om2_hupopdens200), -(om2_hupopdens200)]
dist_lw2_hupopdens200 <- mat2listw(dist_m2_hupopdens200)
moran.test(res.MNE_humanpopdens200, dist_lw2_hupopdens200) # no, I =  -6.398024e-04

res.MNE_humanpopdens100 <- residuals(MNE_humanpopdens100)[,1]
# removing rows from distance matrix because NAs in data
om2_hupopdens100 <-which(is.na(fishy.ne$hupopdens100km), arr.ind=TRUE)
dist_m2_hupopdens100 <- dist_m2[-(om2_hupopdens100), -(om2_hupopdens100)]
dist_lw2_hupopdens100 <- mat2listw(dist_m2_hupopdens100)
moran.test(res.MNE_humanpopdens100, dist_lw2_hupopdens100) # no, I = -6.114707e-04

res.MNE_humanpopdens50 <- residuals(MNE_humanpopdens50)[,1]
# removing rows from distance matrix because NAs in data
om2_hupopdens50 <-which(is.na(fishy.ne$hupopdens50km), arr.ind=TRUE)
dist_m2_hupopdens50 <- dist_m2[-(om2_hupopdens50), -(om2_hupopdens50)]
dist_lw2_hupopdens50 <- mat2listw(dist_m2_hupopdens50)
moran.test(res.MNE_humanpopdens50, dist_lw2_hupopdens50) # no, I = -6.717714e-04

res.MNE_humanpopdens25 <- residuals(MNE_humanpopdens25)[,1]
# removing rows from distance matrix because NAs in data
om2_hupopdens25 <-which(is.na(fishy.ne$hupopdens25km), arr.ind=TRUE)
dist_m2_hupopdens25 <- dist_m2[-(om2_hupopdens25), -(om2_hupopdens25)]
dist_lw2_hupopdens25 <- mat2listw(dist_m2_hupopdens25)
moran.test(res.MNE_humanpopdens25, dist_lw2_hupopdens25) # no, I =  -7.708638e-04 

# Cumulative Human Impacts
res.MNE_huimpacts200 <- residuals(MNE_huimpacts200)[,1]
moran.test(res.MNE_huimpacts200, dist_lw2) # no, I = -7.245401e-04

res.MNE_huimpacts100 <- residuals(MNE_huimpacts100)[,1]
moran.test(res.MNE_huimpacts100, dist_lw2) # no, I =  -7.631160e-04

res.MNE_huimpacts50 <- residuals(MNE_huimpacts50)[,1]
moran.test(res.MNE_huimpacts50, dist_lw2) # no, I = -7.826578e-04

res.MNE_huimpacts25 <- residuals(MNE_huimpacts25)[,1]
moran.test(res.MNE_huimpacts25, dist_lw2) # no, I = -7.688647e-04

# Fishing Pressure
res.MNE_fishing200 <- residuals(MNE_fishing200)[,1]
# removing rows from distance matrix because NAs in data
om2_fishing200 <-which(is.na(fishy.ne$fishing_effort200km), arr.ind=TRUE)
dist_m2_fishing200 <- dist_m2[-(om2_fishing200),-(om2_fishing200)]
dist_lw2_fishing200 <- mat2listw(dist_m2_fishing200)
moran.test(res.MNE_fishing200, dist_lw2_fishing200) # no, I = -7.165941e-04

res.MNE_fishing100 <- residuals(MNE_fishing100)[,1]
# removing rows from distance matrix because NAs in data
om2_fishing100 <-which(is.na(fishy.ne$fishing_effort100km), arr.ind=TRUE)
dist_m2_fishing100 <- dist_m2[-(om2_fishing100),-(om2_fishing100)]
dist_lw2_fishing100 <- mat2listw(dist_m2_fishing100)
moran.test(res.MNE_fishing100, dist_lw2_fishing100) # no, I = -9.006004e-04

res.MNE_fishing50 <- residuals(MNE_fishing50)[,1]
# removing rows from distance matrix because NAs in data
om2_fishing50 <-which(is.na(fishy.ne$fishing_effort50km), arr.ind=TRUE)
dist_m2_fishing50 <- dist_m2[-(om2_fishing50),-(om2_fishing50)]
dist_lw2_fishing50 <- mat2listw(dist_m2_fishing50)
moran.test(res.MNE_fishing50, dist_lw2_fishing50) # no, I = -7.982125e-04 

res.MNE_fishing25 <- residuals(MNE_fishing25)[,1]
# removing rows from distance matrix because NAs in data
om2_fishing25 <-which(is.na(fishy.ne$fishing_effort25km), arr.ind=TRUE)
dist_m2_fishing25 <- dist_m2[-(om2_fishing25),-(om2_fishing25)]
dist_lw2_fishing25 <- mat2listw(dist_m2_fishing25)
moran.test(res.MNE_fishing25, dist_lw2_fishing25) # no, I = -9.803402e-04

# Gene Diversity Models With Priors #

# Human Population Density
res.p.MH_humanpopdens200 <- residuals(p.MH_humanpopdens200)[,1]
moran.test(res.p.MH_humanpopdens200, dist_lw_hupopdens200) # no, I = -3.108130e-03

res.p.MH_humanpopdens100 <- residuals(p.MH_humanpopdens100)[,1]
moran.test(res.p.MH_humanpopdens100, dist_lw_hupopdens100) # no, I = -3.613115e-03

res.p.MH_humanpopdens50 <- residuals(p.MH_humanpopdens50)[,1]
moran.test(res.p.MH_humanpopdens50, dist_lw_hupopdens50) # no, I = -4.827215e-03

res.p.MH_humanpopdens25 <- residuals(p.MH_humanpopdens25)[,1]
moran.test(res.p.MH_humanpopdens25, dist_lw_hupopdens25) # no, I = -4.972154e-03 

# Cumulative Human Impacts
res.p.MH_huimpacts200 <- residuals(p.MH_huimpacts200)[,1]
moran.test(res.p.MH_huimpacts200, dist_lw) # no, I = -4.142218e-03

res.p.MH_huimpacts100 <- residuals(p.MH_huimpacts100)[,1]
moran.test(res.p.MH_huimpacts100, dist_lw) # no, I = -4.389371e-03 

res.p.MH_huimpacts50 <- residuals(p.MH_huimpacts50)[,1]
moran.test(res.p.MH_huimpacts50, dist_lw) # no, I = -4.698327e-03

res.p.MH_huimpacts25 <- residuals(p.MH_huimpacts25)[,1]
moran.test(res.p.MH_huimpacts25, dist_lw) # no, I = -4.747372e-03

# Fishing Pressure
res.p.MH_fishing200 <- residuals(p.MH_fishing200)[,1]
moran.test(res.p.MH_fishing200, dist_lw_fishing200) # no, I = -2.750657e-03

res.p.MH_fishing100 <- residuals(p.MH_fishing100)[,1]
moran.test(res.p.MH_fishing100, dist_lw_fishing100) # no, I = -3.553220e-03

res.p.MH_fishing50 <- residuals(p.MH_fishing50)[,1]
moran.test(res.p.MH_fishing50, dist_lw_fishing50) # no, I = -4.342385e-03

res.p.MH_fishing25 <- residuals(p.MH_fishing25)[,1]
moran.test(res.p.MH_fishing25, dist_lw_fishing25) # no, I = -3.835716e-03

# Allelic Richness Models With Priors #

# Human Population Density
res.p.MAR_humanpopdens200 <- residuals(p.MAR_humanpopdens200)[,1]
moran.test(res.p.MAR_humanpopdens200, dist_lw_hupopdens200) # no, I = -3.152529e-03

res.p.MAR_humanpopdens100 <- residuals(p.MAR_humanpopdens100)[,1]
moran.test(res.p.MAR_humanpopdens100, dist_lw_hupopdens100) # no, I = -3.466101e-03

res.p.MAR_humanpopdens50 <- residuals(p.MAR_humanpopdens50)[,1]
moran.test(res.p.MAR_humanpopdens50, dist_lw_hupopdens50) # no, I = -4.821051e-03

res.p.MAR_humanpopdens25 <- residuals(p.MAR_humanpopdens25)[,1]
moran.test(res.p.MAR_humanpopdens25, dist_lw_hupopdens25) # no, I = -4.850263e-03  

# Cumulative Human Impacts
res.p.MAR_huimpacts200 <- residuals(p.MAR_huimpacts200)[,1]
moran.test(res.p.MAR_huimpacts200, dist_lw) # no, I = -3.632974e-03

res.p.MAR_huimpacts100 <- residuals(p.MAR_huimpacts100)[,1]
moran.test(res.p.MAR_huimpacts100, dist_lw) # no, I = -4.230108e-03

res.p.MAR_huimpacts50 <- residuals(p.MAR_huimpacts50)[,1]
moran.test(res.p.MAR_huimpacts50, dist_lw) # no, I =  -4.578625e-03

res.p.MAR_huimpacts25 <- residuals(p.MAR_huimpacts25)[,1]
moran.test(res.p.MAR_huimpacts25, dist_lw) # no, I = -4.710108e-03  

# Fishing Pressure
res.p.MAR_fishing200 <- residuals(p.MAR_fishing200)[,1]
moran.test(res.p.MAR_fishing200, dist_lw_fishing200) # no, I = -2.795537e-03

res.p.MAR_fishing100 <- residuals(p.MAR_fishing100)[,1]
moran.test(res.p.MAR_fishing100, dist_lw_fishing100) # no, I = -3.777509e-03 

res.p.MAR_fishing50 <- residuals(p.MAR_fishing50)[,1]
moran.test(res.p.MAR_fishing50, dist_lw_fishing50) # no, I =  -4.739269e-03 

res.p.MAR_fishing25 <- residuals(p.MAR_fishing25)[,1]
moran.test(res.p.MAR_fishing25, dist_lw_fishing25) # no, I = -2.627844e-03 

# Fst Models With Priors #

# Human Population Density
res.p.MFST_humanpopdens200 <- residuals(p.MFST_humanpopdens200)[,1]
moran.test(res.p.MFST_humanpopdens200, dist_lw1_hupopdens200) # no, I = -7.270731e-04

res.p.MFST_humanpopdens100 <- residuals(p.MFST_humanpopdens100)[,1]
moran.test(res.p.MFST_humanpopdens100, dist_lw1_hupopdens100) # no, I = -9.086876e-04

res.p.MFST_humanpopdens50 <- residuals(p.MFST_humanpopdens50)[,1]
moran.test(res.p.MFST_humanpopdens50, dist_lw1_hupopdens50) # no, I =  -9.937472e-04 

res.p.MFST_humanpopdens25 <- residuals(p.MFST_humanpopdens25)[,1]
moran.test(res.p.MFST_humanpopdens25, dist_lw1_hupopdens25) # no, I =  -1.067933e-03

# Cumulative Human Impacts
res.p.MFST_huimpacts200 <- residuals(p.MFST_huimpacts200)[,1]
moran.test(res.p.MFST_huimpacts200, dist_lw1) # no, I = -1.043037e-03

res.p.MFST_huimpacts100 <- residuals(p.MFST_huimpacts100)[,1]
moran.test(res.p.MFST_huimpacts100, dist_lw1) # no, I = -1.086263e-03 

res.p.MFST_huimpacts50 <- residuals(p.MFST_huimpacts50)[,1]
moran.test(res.p.MFST_huimpacts50, dist_lw1) # no, I = -1.095671e-03

res.p.MFST_huimpacts25 <- residuals(p.MFST_huimpacts25)[,1]
moran.test(res.p.MFST_huimpacts25, dist_lw1) # no, I = -1.103003e-03

# Fishing Pressure
res.p.MFST_fishing200 <- residuals(p.MFST_fishing200)[,1]
moran.test(res.p.MFST_fishing200, dist_lw1_fishing200) # no, I =  -1.337828e-03 

res.p.MFST_fishing100 <- residuals(p.MFST_fishing100)[,1]
moran.test(res.p.MFST_fishing100, dist_lw1_fishing100) # no, I = -1.464328e-03

res.p.MFST_fishing50 <- residuals(p.MFST_fishing50)[,1]
moran.test(res.p.MFST_fishing50, dist_lw1_fishing50) # no, I = -1.274324e-03

res.p.MFST_fishing25 <- residuals(p.MFST_fishing25)[,1]
moran.test(res.p.MFST_fishing25, dist_lw1_fishing25) # no, I = -8.880658e-04   

# Ne Models With Priors #

# Human Population Density
res.p.MNE_humanpopdens200 <- residuals(p.MNE_humanpopdens200)[,1]
moran.test(res.p.MNE_humanpopdens200, dist_lw2_hupopdens200) # no, I = -6.553328e-04  

res.p.MNE_humanpopdens100 <- residuals(p.MNE_humanpopdens100)[,1]
moran.test(res.p.MNE_humanpopdens100, dist_lw2_hupopdens100) # no, I = -6.110838e-04 

res.p.MNE_humanpopdens50 <- residuals(p.MNE_humanpopdens50)[,1]
moran.test(res.p.MNE_humanpopdens50, dist_lw2_hupopdens50) # no, I = -6.603560e-04 

res.p.MNE_humanpopdens25 <- residuals(p.MNE_humanpopdens25)[,1]
moran.test(res.p.MNE_humanpopdens25, dist_lw2_hupopdens25) # no, I = -7.815503e-04     

# Cumulative Human Impacts
res.p.MNE_huimpacts200 <- residuals(p.MNE_huimpacts200)[,1]
moran.test(res.p.MNE_huimpacts200, dist_lw2) # no, I = -7.188571e-04 

res.p.MNE_huimpacts100 <- residuals(p.MNE_huimpacts100)[,1]
moran.test(res.p.MNE_huimpacts100, dist_lw2) # no, I = -7.554865e-04 

res.p.MNE_huimpacts50 <- residuals(p.MNE_huimpacts50)[,1]
moran.test(res.p.MNE_huimpacts50, dist_lw2) # no, I = -7.724543e-04 

res.p.MNE_huimpacts25 <- residuals(p.MNE_huimpacts25)[,1]
moran.test(res.p.MNE_huimpacts25, dist_lw2) # no, I = -7.624826e-04

# Fishing Pressure
res.p.MNE_fishing200 <- residuals(p.MNE_fishing200)[,1]
moran.test(res.p.MNE_fishing200, dist_lw2_fishing200) # no, I = -7.291080e-04

res.p.MNE_fishing100 <- residuals(p.MNE_fishing100)[,1]
moran.test(res.p.MNE_fishing100, dist_lw2_fishing100) # no, I = -9.035198e-04

res.p.MNE_fishing50 <- residuals(p.MNE_fishing50)[,1]
moran.test(res.p.MNE_fishing50, dist_lw2_fishing50) # no, I = -8.001986e-04

res.p.MNE_fishing25 <- residuals(p.MNE_fishing25)[,1]
moran.test(res.p.MNE_fishing25, dist_lw2_fishing25) # no, I = -1.052656e-03
