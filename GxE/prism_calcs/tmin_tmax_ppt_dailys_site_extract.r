require(reshape2)
require(dplyr)
require(raster)
require(hutilscpp)
require(prism)

#read coordinates
site_coords = read.table("/mnt/home/garnas/ericm/PRISM_analysis/site_coords.USA.txt", header = T)
sites_climate = site_coords


###########################
#Start of loop to loop the files for daily values

prism_set_dl_dir("/mnt/home/garnas/ericm/PRISM_dailies_tmax")
num_files = prism_archive_ls() %>% length

#list for storing data by day
sites_climate_list = list()

for(i in 1:num_files){
    
    print(i)
    new_file<-i#this number corresponds to the row of the file of interest

    #set to tmax
    prism_set_dl_dir("/mnt/home/garnas/ericm/PRISM_dailies_tmax")
    print("tmax")
    
    RS <- pd_stack(prism_archive_ls()[new_file]) ##raster file of data
    proj4string(RS)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
    
    ##convert raster to point data frame
    df <- data.frame(rasterToPoints(RS))
    m.df <- reshape2::melt(df, c("x", "y"))
    names(m.df)[1:2] <- c("lon", "lat") #rename columns
    prism_coords = m.df[,1:2]
    
    #nearest coords#
    
    nrst_coords = match_nrst_haversine(site_coords$lat, site_coords$lon, m.df$lat, m.df$lon)
#    m.df[nrst_coords$pos,]
    
    #Get tmax and date
    sites_climate$tmax = m.df[nrst_coords$pos,colnames(m.df) == "value"]
    coord_day = regmatches(names(RS),regexpr("[0-9]{8,}",names(RS), perl = T)) %>% as.numeric
    sites_climate$day = coord_day
    
    #reset to get tmin
    prism_set_dl_dir("/mnt/home/garnas/ericm/PRISM_dailies_tmin")
    print("tmin")
    
    RS <- pd_stack(prism_archive_ls()[new_file]) ##raster file of data
    proj4string(RS)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
    
    ##convert raster to point data frame
    df <- data.frame(rasterToPoints(RS))
    m.df <- reshape2::melt(df, c("x", "y"))
    names(m.df)[1:2] <- c("lon", "lat") #rename columns
    prism_coords = m.df[,1:2]
    
    #nearest coords#
    
    nrst_coords = match_nrst_haversine(site_coords$lat, site_coords$lon, m.df$lat, m.df$lon)
#    m.df[nrst_coords$pos,]
    
    #Get tmin
    sites_climate$tmin = m.df[nrst_coords$pos,colnames(m.df) == "value"]


    #reset to get ppt
    prism_set_dl_dir("/mnt/home/garnas/ericm/PRISM_dailies_ppt")
    print("ppt")
    
    RS <- pd_stack(prism_archive_ls()[new_file]) ##raster file of data
    proj4string(RS)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

    ##convert raster to point data frame
    df <- data.frame(rasterToPoints(RS))
    m.df <- reshape2::melt(df, c("x", "y"))
    names(m.df)[1:2] <- c("lon", "lat") #rename columns
    prism_coords = m.df[,1:2]

    #nearest coords#

    nrst_coords = match_nrst_haversine(site_coords$lat, site_coords$lon, m.df$lat, m.df$lon)
#    m.df[nrst_coords$pos,]

    #Get ppt
    sites_climate$ppt = m.df[nrst_coords$pos,colnames(m.df) == "value"]

    #Save to list indexed by file number
    sites_climate_list[[i]] = sites_climate

}

#coerce list to df
sites_climate_list.df = do.call(rbind.data.frame, sites_climate_list)

saveRDS(sites_climate_list, "/mnt/home/garnas/ericm/PRISM_analysis/sites_daily_tmin_tmax_ppt_20072018.list.rds")
saveRDS(sites_climate_list.df, "/mnt/home/garnas/ericm/PRISM_analysis/sites_daily_tmin_tmax_ppt_20072018.df.rds")

