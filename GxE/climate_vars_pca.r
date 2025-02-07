#read env data from sample list and unique
library(vegan)

#Join pops to site data

site.info = read.csv("data/sample_metadata/site_info.csv")
site.GDD = read.table("data/sample_metadata/site_climate.GDD.txt", header = T)
site.climate = read.table("data/sample_metadata/sites_climate.txt", header = T)
site.GDD$freezeThaw.annual_mean = site.GDD$freezeThaw.mean_growing + site.GDD$freezeThaw.mean_nongrowing

site_metadata = left_join(site.GDD, site.info %>% select(Site, lat, lon, duration_infection), by = "Site") %>%
    left_join(., site.climate %>% select(Site, tmin, tmax, ppt, MAT, lat, lon, state.name), by = c("Site", "lat", "lon") )


sample_metadata.site_info = left_join(sample_metadata, site_metadata, by = "state.name")
colnames(sample_metadata.site_info)

sample_metadata.site_info.uniq = sample_metadata.site_info %>% select(-sample) %>% distinct()

#scale the vars

sample_metadata.scaled = apply(
    sample_metadata.site_info.uniq %>% 
        select(HDD4.mean_growing, HDD4.mean_nongrowing, tmin, ppt, freezeThaw.annual_mean), 
    2, 
    scale
)

#run princomp

climate.pca = capscale(sample_metadata.scaled ~ 1, distance = "euclidean")

str(climate.pca)
plot(climate.pca$CA$eig/sum(climate.pca$CA$eig) )
#First three axes contain majority of var
biplot(climate.pca)
biplot(climate.pca, choices = c(1,3))
biplot(climate.pca, choices = c(2,3))

pdf("figures/climate_vars_pca.pdf")
plot(climate.pca$CA$eig/sum(climate.pca$CA$eig) )
biplot(climate.pca)
biplot(climate.pca, choices = c(1,3))
biplot(climate.pca, choices = c(2,3))
dev.off()
