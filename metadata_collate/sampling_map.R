library(ggplot2)
library(sf)
library(stringr)
library(dplyr)
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
    filter(PRVSTTNAME %in% c(scale_dat$PRVSTTNAME, "Kentucky", "Illinois", "Indiana", "Georgia", "Alabama", "South Carolina", "Iowa", "Mississippi", "Arkansas", "Missouri"))
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

meta_n = rbind(Nf.meta_n, Nd.meta_n)
meta_n %>% print(n = Inf)
#we're going to move the MI and MI.UP a bit manuualy so both sets are visible (they have the same n)
#Nf went left by .1 and Nd went right by .1
meta_n[4,3] = -84.8
meta_n[5,3] = -87.2
meta_n[26,3] = -84.6
meta_n[27,3] = -87.0
#same deal for the Odell Nf/Nd (NB.YO)
meta_n[8,3] = -66.8
meta_n[30,3] = -66.6
#also for the Gatineau Nf/Nd (QC.OUG), BUT they are also overplotting on Chelsea (QC.OUG) which has three Nf
# We are going to need to rerun the within between tests too, because we had the QC.OU sites as the same (but they are not)
#Nf went left by .1 and Nd went right by .1
meta_n[22,3] = -75.8
meta_n[41,3] = -76.0
#also move down slightly
meta_n[22,2] = 45.19
meta_n[41,2] = 45.21




## there is a lot of overlap in the NH site plotting. Let's move CW a bit to the right, move JG a bit south and right
meta_n %>% print(n = Inf)
# NH.JG    43.1 -70.9
meta_n[32,2] = 43.1
meta_n[32,3] = -70.725
# NH.CW    43.1 -71.0
meta_n[12,3] = -70.5
#move wf left and down
#NH.WF    43.2 -70.9     1 Nd   
meta_n[34,3] = -71.15
meta_n[34,2] = 43.1
#NH.WF    43.2 -70.9     5 Nf   
meta_n[14,3] = -71.3
meta_n[14,2] = 43.1
# 31 NH.CCM   43.5 -71.1     6 Nd   
# 33 NH.STR   43.2 -71.1     1 Nd  
# 11 NH.CCM   43.5 -71.1     1 Nf 
# let's move STR left and CCM Nf right
#meta_n[11,3] = -71.1
#meta_n[11,2] = 43.5
meta_n[33,3] = -71.475
meta_n[33,2] = 43.2
# the Nf point is plotting below the Nd.
# might just need to annotate on top of it

p1 = ggplot(all.counties.scale) +
    geom_sf(aes(fill = SCALEYR)) +
    scale_fill_gradient2(low = "#2166ac", high = "#b2182b", mid = "#d1e5f0", midpoint = 1950) +
    coord_sf(default_crs = sf::st_crs(4326), xlim = c(-89.25,-61), ylim = c(34, 50.25)) +
    geom_point(
        data = meta_n, 
        aes(
            x = lon, 
            y = lat, 
            color = factor(spp, levels = c("Nf", "Nd"), labels = c("N. faginata", "N. ditissima")), 
            size = n
        )
    ) +
    geom_point(data = meta_n %>% filter(spp != "Nf"),  
        aes(
            x = lon, 
            y = lat, 
            size = n
        ),
        shape = 1,
        colour = "black"
    ) +
    annotate(geom = "point", x = -71.1, y = 43.5, size = 1.9, color = "black") +
    scale_size_continuous(breaks = c(1,3,6,9,12), range = c(2, 7)) +
    labs(color = "Species", size = "Individuals (n)", fill = "Year scale detected") +
    scale_color_manual(values = c("black", "white")) +
    my_gg_theme.def_size +
    guides(
        color = guide_legend(
            order = 1,
            override.aes = list(
                fill = c("black", "white"), 
                shape = 21, 
                color = "black", 
                size = 3
            )
        ),
        size = guide_legend(order = 2),
        fill = guide_colorbar(order = 3)
    ) +
    theme(
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        legend.key = element_blank()
    )
p1

pdf("figures/pop_gen/sample_n_map.pdf", width = 8, height = 7)
p1
dev.off()

