library(geodata)
library(dplyr)
library(dismo)

# we first donwload the required monthly data and filter to coords
# we then use dismo::biovars to calculate the bioclim vars


# metadata for filtering to 
Nf.meta = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
Nd.meta = read.csv("data/sample_metadata/Nd_filtered.lat_lon_dur_inf.csv")

Nf.meta_n = Nf.meta %>% 
    dplyr::select(state, lat, lon) %>% 
    group_by(state, lat, lon) %>% 
    summarize(n = n())
Nf.meta_n$spp = "Nf"
Nd.meta_n = Nd.meta %>% 
    dplyr::select(state, lat, lon) %>% 
    group_by(state, lat, lon) %>% 
    summarize(n = n())
Nd.meta_n$spp = "Nd"

# the coords for Nf
coords = data.frame(lon = Nf.meta_n$lon, lat = Nf.meta_n$lat)
## Convert coords data frame to vector object
points <- vect(coords,
               geom=c("lon", "lat"),
               crs = "EPSG:4326"
    )

# download tmin to test
tmin = worldclim_tile(var = "tmin", 
    lon = Nf.meta$lon[1], 
    lat = Nf.meta$lat[1], 
    path = "data/sample_metadata/world_clim", 
    version="2.1", 
    res = 0.5
)
tmin 

# the tile this grabs runs from -60 to -90 lon and from 30 to 60 lat
range(coords$lon)
range(coords$lat)
#just makes it!
#

tmin.values = extract(
    tmin,
    points,
    ID = F,
    xy = T
)
tmin.values
class(tmin.values)
#the  result is a df of one row per coord with 1 column per month
# the input to bio vars is vector of 12 months precip, tmn, tmx as below so 
# can input rows consecutively


tmax = worldclim_tile(var = "tmax", 
    lon = Nf.meta$lon[1], 
    lat = Nf.meta$lat[1], 
    path = "data/sample_metadata/world_clim", 
    version="2.1", 
    res = 0.5
)

prec = worldclim_tile(var = "prec", 
    lon = Nf.meta$lon[1], 
    lat = Nf.meta$lat[1], 
    path = "data/sample_metadata/world_clim", 
    version="2.1", 
    res = 0.5
)

# must use vector as input. The func stupidly gives different output based on 
# class (and does not give clear explanation as to what)
#biovars(as.numeric(prec.values[1,1:12]), as.numeric(tmin.values[1,1:12]), as.numeric(tmax.values[1,1:12])) %>% data.frame()





# we do the point set and extraction line by line to 
# 1. maintain order and 
# 2. make this more modular 

biovars_list = list()
for(i in 1:length(Nf.meta_n$state)){
    
    points = vect(
        Nf.meta_n[i,c("lon", "lat")],
        geom=c("lon", "lat"),
        crs = "EPSG:4326"
    )
    
    tmin.values = extract(
        tmin,
        points,
        ID = F,
        xy = F
    )
    tmax.values = extract(
        tmax,
        points,
        ID = F,
        xy = F
    )
    prec.values = extract(
        prec,
        points,
        ID = F,
        xy = F
    )
    biovars_list[[ Nf.meta_n$state[i] ]] = data.frame(
        state = Nf.meta_n$state[i],
        lat = Nf.meta_n$lat[i],
        lon = Nf.meta_n$lon[i],
        biovars(
            as.numeric(prec.values), 
            as.numeric(tmin.values), 
            as.numeric(tmax.values)
        )
    )
}

Nf.biovars = bind_rows(biovars_list)


biovars_list = list()
for(i in 1:length(Nd.meta_n$state)){
    
    points = vect(
        Nd.meta_n[i,c("lon", "lat")],
        geom=c("lon", "lat"),
        crs = "EPSG:4326"
    )
    
    tmin.values = extract(
        tmin,
        points,
        ID = F,
        xy = F
    )
    tmax.values = extract(
        tmax,
        points,
        ID = F,
        xy = F
    )
    prec.values = extract(
        prec,
        points,
        ID = F,
        xy = F
    )
    biovars_list[[ Nd.meta_n$state[i] ]] = data.frame(
        state = Nd.meta_n$state[i],
        lat = Nd.meta_n$lat[i],
        lon = Nd.meta_n$lon[i],
        biovars(
            as.numeric(prec.values), 
            as.numeric(tmin.values), 
            as.numeric(tmax.values)
        )
    )
}

Nd.biovars = bind_rows(biovars_list)
# we tried close points for NS.SE and it does not exist

write.csv(Nf.biovars, "data/sample_metadata/Nf.site_bioclim.csv", row.names = F, quote = F)
write.csv(Nd.biovars, "data/sample_metadata/Nd.site_bioclim.csv", row.names = F, quote = F)
