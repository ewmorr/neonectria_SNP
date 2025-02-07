library(devtools)
library(reshape2)
library(dplyr)
library(raster)

library(hutilscpp)
#install_github(repo = "ropensci/prism")
library(prism)

prism_set_dl_dir("/Users/ericmorrison/PRISM_data")
#options(prism.path = "/Users/ericmorrison/PRISM_data")
#get_prism_normals(type = 'ppt', resolution = '4km', annual = T, keepZip = TRUE)
#get_prism_normals(type = 'tmean', resolution = '4km', annual = T, keepZip = TRUE)

prism_set_dl_dir("/Users/ericmorrison/PRISM_data_monthly")
#get_prism_normals(type = 'tmax', resolution = '4km', annual = F, keepZip = TRUE, mon = c(6,7,8,9))
#get_prism_normals(type = 'tmin', resolution = '4km', annual = F, keepZip = TRUE, mon = c(1,2,3,12))

# metadata
site_coords = read.csv("data/sample_metadata/site_info.csv")
sites_climate = site_coords %>% filter(country == "USA")
write.table(sites_climate %>% dplyr::select(Site, lat, lon, state), "data/sample_metadata/site_coords.USA.txt", quote = F, sep = "\t", row.names = F)

#####
#ppt#
prism_set_dl_dir("/Users/ericmorrison/PRISM_data")

new_file<-1#this number corresponds to the row of the file of interest
RS <- pd_stack(prism_archive_ls()[new_file]) ##raster file of data
proj4string(RS)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

####################################
##convert raster to point data frame
df <- data.frame(rasterToPoints(RS))
m.df <- reshape2::melt(df, c("x", "y"))
names(m.df)[1:2] <- c("lon", "lat") #rename columns
prism_coords = m.df[,1:2]
################
#nearest coords#
nrst_coords = match_nrst_haversine(sites_climate$lat, sites_climate$lon, m.df$lat, m.df$lon)
m.df[nrst_coords$pos,]
####################################

#get ppt of nearest
sites_climate$ppt = m.df[nrst_coords$pos,colnames(m.df) == "value"]

#####
#MAT#
new_file<-3#this number corresponds to the row of the file of interest
RS <- pd_stack(prism_archive_ls()[new_file]) ##raster file of data
proj4string(RS)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

##convert raster to point data frame
df <- data.frame(rasterToPoints(RS))
m.df <- reshape2::melt(df, c("x", "y"))
names(m.df)[1:2] <- c("lon", "lat") #rename columns


sites_climate$MAT = m.df[nrst_coords$pos,colnames(m.df) == "value"]


#####
#tmax#
prism_set_dl_dir("/Users/ericmorrison/PRISM_data_monthly")

temp.df = data.frame(Site = sites_climate$Site)

for(i in 1:4){
    
    new_file <- i
    RS <- pd_stack(prism_archive_ls()[new_file]) ##raster file of data
    proj4string(RS)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

    ##convert raster to point data frame
    df <- data.frame(rasterToPoints(RS))
    m.df <- reshape2::melt(df, c("x", "y"))
    names(m.df)[1:2] <- c("lon", "lat") #rename columns
    
    temp.tmax = data.frame(m.df[nrst_coords$pos,colnames(m.df) == "value"])
    colnames(temp.tmax) = paste("month", i, sep = ".")
    temp.df = cbind(temp.df, temp.tmax)
}
temp.df
#July is the hottest month at each site, but barely. We will get max by row
apply(temp.df[2:5], 1, max)

sites_climate$tmax = apply(temp.df[2:5], 1, max)


#####
#tmin#

temp.df = data.frame(Site = sites_climate$Site)

for(i in 5:8){
    
    new_file <- i
    RS <- pd_stack(prism_archive_ls()[new_file]) ##raster file of data
    proj4string(RS)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
    
    ##convert raster to point data frame
    df <- data.frame(rasterToPoints(RS))
    m.df <- reshape2::melt(df, c("x", "y"))
    names(m.df)[1:2] <- c("lon", "lat") #rename columns
    
    temp.tmin = data.frame(m.df[nrst_coords$pos,colnames(m.df) == "value"])
    colnames(temp.tmin) = paste("month", i, sep = ".")
    temp.df = cbind(temp.df, temp.tmin)
}
temp.df
#Jan is the coldest month in all but 2 sites. We will get min by row
apply(temp.df[2:5], 1, min)

sites_climate$tmin = apply(temp.df[2:5], 1, min)

#######
#write#

write.table(sites_climate, file = "data/sample_metadata/sites_climate.txt", quote = F, sep = "\t", row.names = F)
