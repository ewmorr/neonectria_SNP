library(ggplot2)
library(sf)
library(stringr)
library(dplyr)
library(tidyr)
library(scatterpie)
source("library/ggplot_theme.txt")

#map data from here https://exploratory.io/map
#https://blog.exploratory.io/making-maps-for-canadas-provisions-and-census-divisions-in-r-c189b88ccd8a

#counties = FROM_GeoJson(url_file_string = "data/sample_metadata/county/county.geojson")
counties = st_read("data/sample_metadata/county/county.geojson")
districts = st_read("data/sample_metadata/canada_divisions/canada_divisions.geojson")
head(counties)
head(districts)

districts$COUNTRY = "Canada"
counties$COUNTRY = "USA"
#remove french prov names
districts$PRNAME = sub(" / .*","",districts$PRNAME)
#split county state name
county_state = str_split_fixed(counties$COUNTY_STATE_NAME, ", ", n = 2)
counties$county = county_state[,1]
counties$state = county_state[,2]

districts.cut = districts[,c("COUNTRY", "PRNAME", "CDNAME", "geometry")]
counties.cut = counties[,c("COUNTRY", "state", "county", "geometry")]
#name to the same as the Cale durInf file
colnames(districts.cut) = c("COUNTRY", "PRVSTTNAME", "CONAME","geometry")
colnames(counties.cut) = c("COUNTRY", "PRVSTTNAME", "CONAME","geometry")
all.counties = rbind(districts.cut, counties.cut)

########################
# scale infestation data
scale_dat = read.csv("data/sample_metadata/Cale_Morin-BeechScaleDatesCanadaUS.csv")
scale_dat = scale_dat[,1:4]
head(scale_dat)
all.counties.scale_slice = all.counties %>% 
    filter(PRVSTTNAME %in% c(scale_dat$PRVSTTNAME, "Kentucky", "Illinois", "Indiana", "Georgia", "Alabama", "South Carolina", "Iowa", "Mississippi", "Arkansas", "Missouri", "Delaware"))
nrow(all.counties.scale_slice)
nrow(scale_dat)
#check to see if we pick up all the counties
all.counties %>% filter(PRVSTTNAME %in% scale_dat$PRVSTTNAME & CONAME %in% scale_dat$CONAME) %>% nrow
# it.s more than the rows of scale_dat... that's odd
# 
#We did not get plotting for a lot of Quebec. Looks like French accent characters didn't convert correctly. 
#There are also more Quebeclisted in the scale data than in the subdiv data (which explains the
# discrepancy in nrow above
districts.cut %>% filter(PRVSTTNAME == "Quebec")
scale_dat %>% filter(PRVSTTNAME == "Quebec")
# note that the subdiv data I have are much higher resolution than what is plotted in the Cale paper
# There are 98 Quebec features in the subdiv data and 14 in the Cale data
districts.cut %>% filter(PRVSTTNAME == "Quebec") %>% pull(CONAME)
scale_dat %>% filter(PRVSTTNAME == "Quebec") %>% pull(CONAME)
#We checked the province level data from the same source but it doesn't help.
#Looks like we will have to do this manually
subdiv_names.quebec = districts.cut %>% filter(PRVSTTNAME == "Quebec") %>% pull(CONAME)
subdivScale_names.quebec = scale_dat %>% filter(PRVSTTNAME == "Quebec") %>% pull(CONAME)

#write.table(sort(subdiv_names.quebec), "data/sample_metadata/quebec_bad_names.txt", row.names = F, col.names = F, quote = F)
#write.table(sort(subdivScale_names.quebec), "data/sample_metadata/quebec_scale_good_names.txt", row.names = F, col.names = F, quote = F)

#manually mapped the subdiv names to region names based on here
# https://en.wikipedia.org/wiki/List_of_regional_county_municipalities_and_equivalent_territories_in_Quebec
#Nord-du-quebec_removed (it's not in the scale data and we don't map it
quebec_name_map = read.csv("data/sample_metadata/quebec_bad_names_mapping.txt", header = F, sep = "\t")
head(scale_dat)
scale_dat.qc = scale_dat %>% filter(PRVSTTNAME == "Quebec")
colnames(scale_dat.qc)[3] = "region"
colnames(quebec_name_map) = c("CONAME", "region")
scale_dat.qc_new = left_join(quebec_name_map, scale_dat.qc) %>% select(COUNTRY, PRVSTTNAME, CONAME, SCALEYR)

scale_dat.new = rbind(scale_dat %>% filter(PRVSTTNAME != "Quebec"), scale_dat.qc_new)
nrow(scale_dat.new)
#join the data
all.counties.scale = left_join(all.counties.scale_slice, scale_dat.new, by = c("COUNTRY", "PRVSTTNAME", "CONAME"))
all.counties.scale %>% filter(PRVSTTNAME == "Quebec") 
# mapping worked for the accent charcters, nice

#read sample metadata for plotting sample numbers
Nf.meta = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
Nd.meta = read.csv("data/sample_metadata/Nd_filtered.lat_lon_dur_inf.csv")
Nf.meta %>% filter(state == "VA")
all.counties.scale %>% filter(PRVSTTNAME == "Virginia" & CONAME == "Rockingham")
#the cale dat has 2014 and we have 2013. We will change since we have NP data
all.counties.scale[all.counties.scale$PRVSTTNAME == "Virginia" & all.counties.scale$CONAME == "Rockingham", "SCALEYR"] = 2013
all.counties.scale %>% filter(PRVSTTNAME == "Virginia" & CONAME == "Rockingham")


Nf.meta %>% 
    group_by(state, lat, lon, duration_infection, Site, collection_period) %>% 
    summarize(n = n()) %>%
    print(n = Inf)
Nf.meta %>% 
    group_by(state, lat, lon, duration_infection, Site, collection_period) %>% 
    summarize(n = n()) %>%
    pull(n) %>% sum()
Nd.meta %>% 
    group_by(state, lat, lon, duration_infection, Site, collection_period) %>% 
    summarize(n = n()) %>%
    print(n = Inf)
Nd.meta %>% 
    group_by(state, lat, lon, duration_infection,collection_period) %>% 
    summarize(n = n()) %>%
    print(n = Inf)
Nd.meta %>% 
    group_by(state, lat, lon, duration_infection, Site, collection_period) %>% 
    summarize(n = n()) %>%
    pull(n) %>% sum()

Nf.meta_n = Nf.meta %>% 
    select(state, lat, lon) %>% 
    group_by(state, lat, lon) %>% 
    summarize(n = n())
Nf.meta_n$spp = "Nf"
Nd.meta_n = Nd.meta %>% 
    select(state, lat, lon) %>% 
    group_by(state, lat, lon) %>% 
    summarize(n = n())
Nd.meta_n$spp = "Nd"

#################################
#admixture dat Nf K = 2, Nd K = 5
#metadata and ancestry 
Nf.fam_info = read.table("data/Nf/final_tables/rm_dups/FINAL_snp.admixture.nosex", header = F)
colnames(Nf.fam_info)[1] = "Sequence_label"
Nf.sample_metadata = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
Nf.sample_metadata = left_join(data.frame(Sequence_label = Nf.fam_info[,1]), Nf.sample_metadata)

Nd.fam_info = read.table("data/Nd/final_tables/rm_dups/FINAL_snp.admixture.nosex", header = F)
colnames(Nd.fam_info)[1] = "Sequence_label"
Nd.sample_metadata = read.csv("data/sample_metadata/Nd_filtered.lat_lon_dur_inf.csv")
Nd.sample_metadata = left_join(data.frame(Sequence_label = Nd.fam_info[,1]), Nd.sample_metadata)


Nf.K_2 = read.table("data/Nf/admixture/FINAL_snp.admixture.2.Q", header = F)
colnames(Nf.K_2) = paste0("Q", 1:ncol(Nf.K_2))
Nf.K_2$Sequence_label = Nf.fam_info$Sequence_label
Nf.K_2$K = 2
Nf.K_2 = Nf.K_2  %>% 
    pivot_longer(cols = starts_with("Q"), names_to = "ancestor", values_to = "Q")

Nd.K_5 = read.table("data/Nd/admixture/FINAL_snp.admixture.5.Q", header = F)
colnames(Nd.K_5) = paste0("Q", 1:ncol(Nd.K_5))
Nd.K_5$Sequence_label = Nd.fam_info$Sequence_label
Nd.K_5$K = 5
Nd.K_5 = Nd.K_5  %>% 
    pivot_longer(cols = starts_with("Q"), names_to = "ancestor", values_to = "Q")

Nf.K_2.meta = left_join(Nf.K_2, Nf.sample_metadata %>% select(Sequence_label, state), by = "Sequence_label") 
nrow(Nf.K_2.meta)
head(Nf.K_2.meta)
Nf.K_2.meta %>%
    ungroup() %>%
    group_by(ancestor, state) %>%
    summarise(mean_Q = mean(Q)) -> Nf.K_2.site_mean

Nd.K_5.meta = left_join(Nd.K_5, Nd.sample_metadata %>% select(Sequence_label, state), by = "Sequence_label") 
nrow(Nd.K_5.meta)
head(Nd.K_5.meta)
Nd.K_5.meta %>%
    ungroup() %>%
    group_by(ancestor, state) %>%
    summarise(mean_Q = mean(Q)) -> Nd.K_5.site_mean

Nf.K_2.wide = Nf.K_2.site_mean %>% pivot_wider(names_from = ancestor, values_from = mean_Q)
Nd.K_5.wide = Nd.K_5.site_mean %>% pivot_wider(names_from = ancestor, values_from = mean_Q)
#################################

Nf.K_2.meta.map = left_join(Nf.K_2.wide, Nf.meta_n)
Nd.K_5.meta.map = left_join(Nd.K_5.wide, Nd.meta_n)
nrow(Nf.K_2.meta.map)
nrow(Nd.K_5.meta.map)

Nf.K_2.meta.map$nmax = ifelse(Nf.K_2.meta.map$n > 4, 4, Nf.K_2.meta.map$n)
Nf.K_2.meta.map$nmax_scaled = log10((Nf.K_2.meta.map$nmax/max(Nf.K_2.meta.map$nmax))+1)

log10(Nf.K_2.meta.map$nmax_scaled+1)

Nf.K_2.meta.map %>% print(n = Inf)
# push NH.CW left and NH. wf right
# 12 NH.CW   0.660   0.340  43.1 -71.0    11 Nf   
# 14 NH.WF   0.0371  0.963  43.2 -70.9     5 Nf   
Nf.K_2.meta.map[12,"lon"] = -71.2
Nf.K_2.meta.map[14,"lon"] = -70.7
Nf.K_2.meta.map[12,"lat"] = 42.9
Nf.K_2.meta.map[14,"lat"] = 42.9

#21 QC.OUC  0.0430  0.957  45.5 -75.8     3 Nf        3      0.243 
#22 QC.OUG  0.0434  0.957  45.5 -75.9     1 Nf        1      0.0969
Nf.K_2.meta.map[21,"lon"] = -75.5
Nf.K_2.meta.map[22,"lon"] = -76.2

#6 NB.LUD  0.0228  0.977  46.5 -66.4     1 Nf        1      0.0969
#7 NB.RB   0.0624  0.938  46.8 -66.4     1 Nf        1      0.0969
#8 NB.YO   0.00427 0.996  46.0 -66.7     1 Nf        1      0.0969
Nf.K_2.meta.map[8,"lat"] = 45.8
Nf.K_2.meta.map[6,"lat"] = 46.2
Nf.K_2.meta.map[6,"lon"] = -66.2

#this is default
default_crs = sf::st_crs(4326)
# this is NAD83
default_crs = sf::st_crs(4269)

st_crs(all.counties.scale)

# because scatterpie has a bug that causes it to plot ovals with geom_sf
# we extract the data to a df and plot with geom_polygon
all.counties.coords = data.frame(sf::st_coordinates(all.counties.scale))
all.counties.coords$id = paste0(all.counties.coords[,3],all.counties.coords[,4],all.counties.coords[,5])
colnames(all.counties.coords)[1:2] = c("x", "y")

p10 = ggplot(all.counties.coords) +
    geom_polygon(aes(x = x, y = y, group = id), fill = "#636363") +
    coord_fixed(xlim = c(-87.5,-61), ylim = c(34.5, 49), ratio = 1) +
    geom_scatterpie2(data = Nf.K_2.meta.map, 
        aes(x = lon, y = lat), 
        cols = paste0("Q", 1:2)#,
#        pie_scale = 1
    ) +
    scale_fill_brewer(palette = "Set1", guide = "none") +
    my_gg_theme.def_size +
    labs(title = "f") +
    theme(
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        legend.key = element_blank(),
        plot.title = element_text(hjust = -0.0175, vjust = -1.5, size = 20),
        plot.margin = margin(t = -10, r = 5.5, b = 5.5, l = 5.5)
    )
p10
    
saveRDS(p10, "figures/pop_gen/admixture/Nf.avg_map.rds")

Nd.K_5.meta.map %>% print(n = Inf)
#7 NH.CCM 0.0000102 0.0362   0.741   0.0447  0.178    43.5 -71.1     6 Nd   
#8 NH.JG  0.00001   0.00001  0.00001 1.00    0.00001  43.1 -70.9     1 Nd   
#9 NH.STR 0.00001   0.00001  1.00    0.00001 0.00001  43.2 -71.1     1 Nd   
#10 NH.WF  0.00759   0.223    0.410   0.307   0.0521   43.2 -70.9     1 Nd   
Nd.K_5.meta.map[8,"lat"] = 42.8
Nd.K_5.meta.map[8,"lon"] = -70.5
Nd.K_5.meta.map[9,"lat"] = 43.3
Nd.K_5.meta.map[9,"lon"] = -71.7
Nd.K_5.meta.map[7,"lat"] = 44.1
Nd.K_5.meta.map[10,"lat"] = 43.4

 
#5 NB.NEW 0.00001   0.00001  0.00001 1.00    0.00001  47.2 -65.6     1 Nd   
#6 NB.YO  1.00      0.000013 0.00001 0.00001 0.00001  46.0 -66.7     1 Nd   

#4 NB.AL  0.00001   0.00001  0.00001 1.00    0.00001  45.5 -65.2     1 Nd  
#11 NS.CU  0.00001   0.00001  0.00001 1.00    0.00001  45.4 -64.9     2 Nd   


#12 NS.GU  0.000016  0.00001  0.00001 0.00001 1.00     45.2 -62.0     1 Nd   

Nd.K_5.meta.map[4,"lon"] = -65.4
Nd.K_5.meta.map[11,"lon"] = -64.7


p11 = ggplot(all.counties.coords) +
    geom_polygon(aes(x = x, y = y, group = id), fill = "#636363") +
    coord_fixed(xlim = c(-87.5,-61), ylim = c(34.5, 49), ratio = 1) +
    geom_scatterpie2(data = Nd.K_5.meta.map, 
        aes(x = lon, y = lat), 
        cols = paste0("Q", 1:5)#,
#        pie_scale = 1
    ) +
    scale_fill_brewer(palette = "Paired", guide = "none") +
    my_gg_theme.def_size +
    labs(title = "g") +
    theme(
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        legend.key = element_blank(),
        plot.title = element_text(hjust = -0.02, vjust = -1.5, size = 20),
        plot.margin = margin(t = -10, r = 5.5, b = 5.5, l = 5.5)
    )
p11
    

saveRDS(p11, "figures/pop_gen/admixture/Nd.avg_map.rds")

