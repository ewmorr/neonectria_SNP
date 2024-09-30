library(dplyr)
library(ggplot2)
library(gridExtra)
library(vegan)
library(geosphere)
source("library/ggplot_theme.txt")

sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
sample_metadata.Nd = read.csv("data/sample_metadata/Nd_filtered.lat_lon_dur_inf.csv")

sample_metadata = rbind(sample_metadata.Nf %>% filter(collection_period == "modern"), 
    sample_metadata.Nd %>% filter(collection_period == "modern")
) %>% select(duration_infection, lat, lon, state) %>%
    distinct

Dgeo = distm(x = sample_metadata[,c("lon", "lat")], fun = distVincentyEllipsoid)
#rownames(Dgeo) = sample_metadata$state
#colnames(Dgeo) = sample_metadata$state
Ddur = dist(x = sample_metadata$duration_infection) %>% as.matrix
#rownames(Ddur) = sample_metadata$state
#colnames(Ddur) = sample_metadata$state

#we are not naming by state because there are two QC.OU sites with different coords
#this is because there is one for Nf and one for Nd

mantel(Ddur, Dgeo)
#Mantel statistic r:     0.7375 
#Significance: 0.001 

Dgeo.long = reshape2::melt(Dgeo)
Dgeo.long %>% filter(Var1 == "NY.N" & Var2 == "NY.S")
Ddur.long = reshape2::melt(Ddur)
colnames(Dgeo.long)[3] = "km"
Dgeo.long$km = Dgeo.long$km/1000
colnames(Ddur.long)[3] = "durDif"
#Dgeo.long$joinby = paste0(Dgeo.long$Var1, Dgeo.long$Var2)
#Ddur.long$joinby = paste0(Ddur.long$Var1, Ddur.long$Var2)

full_dat = left_join(Dgeo.long, Ddur.long, by = c("Var1", "Var2"))

p1 = ggplot(full_dat, aes(x = km, y = durDif)) +
    geom_point(alpha = 0.25, shape = 1) +
    geom_smooth(method = "lm", linetype = 2, color = "black") +
    labs(x = "Geographic distance (km)", y = "Difference in infestation duration (years)") +
    annotate(
        geom = "text", 
        label = expression(paste("Mantel r = 0.74,", italic("P"), " = 0.001")),
        x = 1750,
        y = 2.5
    ) +
    my_gg_theme.def_size 
p1

pdf("figures/pop_gen/IBD/geo_v_infestation_duration.pdf", width = 6.5, height = 3.5)
p1
dev.off()

