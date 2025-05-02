#read env data from sample list and unique
library(vegan)
library(dplyr)
library(ggplot2)
library(GGally)
source("library/ggplot_theme.txt")

#Join pops to site data
site.info = read.csv("data/sample_metadata/site_info.csv")
site.info = site.info %>% filter(country == "USA")
site.GDD = read.table("data/PRISM_dailies/site_climate.GDD.txt", header = T)
site.climate = read.csv("data/sample_metadata/sites_climate.csv", header = T)
site.GDD$freezeThaw.annual_mean = site.GDD$freezeThaw.mean_growing + site.GDD$freezeThaw.mean_nongrowing

site_metadata = left_join(site.GDD, site.info %>% select(Site, lat, lon, duration_infection), by = "Site") %>%
    left_join(., site.climate %>% select(Site, tmin, tmax, ppt, MAT, lat, lon, state), by = c("Site", "lat", "lon") ) %>%
    select(-Site) %>%
    unique()
site_metadata = site_metadata[-14,] #remove extra nh.ccm


#sample_metadata.site_info = left_join(sample_metadata, site_metadata, by = "state")
#colnames(sample_metadata.site_info)
#sample_metadata.site_info.uniq = sample_metadata.site_info %>% select(-sample) %>% distinct()

#scale the vars

site_metadata.scaled = apply(
    site_metadata %>% 
        select(HDD4.mean_growing, HDD4.mean_nongrowing, tmin, tmax, MAT, ppt, freezeThaw.annual_mean, duration_infection), 
    2, 
    scale
) %>% as.data.frame
rownames(site_metadata.scaled) = site_metadata$state
head(site_metadata)

#run princomp
climate.pca = capscale(site_metadata.scaled ~ 1, distance = "euclidean")

str(climate.pca)
plot(climate.pca$CA$eig/sum(climate.pca$CA$eig) )
climate.pca$CA$eig/sum(climate.pca$CA$eig)
#        MDS1         MDS2         MDS3         MDS4         MDS5         MDS6         MDS7         MDS8 
#0.5398403578 0.2074280322 0.1248073735 0.1105816920 0.0105971982 0.0046546423 0.0016979594 0.0003927446 
sum((climate.pca$CA$eig/sum(climate.pca$CA$eig))[1:2])
#First two axes contain  0.7472684
sum((climate.pca$CA$eig/sum(climate.pca$CA$eig))[1:3])
#First three axes contain 0.8720758 of var
biplot(climate.pca)
biplot(climate.pca, choices = c(1,3))
biplot(climate.pca, choices = c(2,3))

#
pca.sites = data.frame(scores(climate.pca, choices = 1:3)$sites)
pca.spp = data.frame(scores(climate.pca, choices = 1:3)$species)
pca.sites$state = rownames(pca.sites)
pca.spp$var = rownames(pca.spp)
pca.spp$var = c(
    "GDD (growing)",
    "GDD (nongrowing)",
    "Tmin",
    "Tmax",
    "MAT",
    "PPT",
    "freeze-thaw",
    "Infestation duration"
)

site_cols.df = read.csv("data/sample_metadata/Nf_site_colors.csv")
pca.sites = left_join(pca.sites, site_cols.df)
site_cols = pca.sites$colors
names(site_cols) = pca.sites$state

p1 = ggplot() +
    geom_hline(yintercept = 0, lty = 2) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_segment(data = pca.spp, 
        aes(x = rep(0, nrow(pca.spp)), y = rep(0, nrow(pca.spp)), xend = MDS1, yend = MDS2),
        arrow = arrow(length = unit(0.1, "cm"))
    ) +
    geom_text(data = pca.spp, aes(x = MDS1*1.1, y = MDS2*1.1, label = var)) +
    geom_point(data = pca.sites, aes(x = MDS1, y = MDS2, color = state), size = 3) +
    scale_color_manual(values = site_cols) +
    labs(x = "PC1 (54%)", y = "PC2 (21%)", color = "Site") +
    my_gg_theme.def_size +
    theme(
        legend.title = element_text(size = 12)
    )

p2 = ggplot() +
    geom_hline(yintercept = 0, lty = 2) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_segment(data = pca.spp, 
        aes(x = rep(0, nrow(pca.spp)), y = rep(0, nrow(pca.spp)), xend = MDS1, yend = MDS3),
        arrow = arrow(length = unit(0.1, "cm"))
    ) +
    geom_text(data = pca.spp, aes(x = MDS1*1.1, y = MDS3*1.1, label = var)) +
    geom_point(data = pca.sites, aes(x = MDS1, y = MDS3, color = state), size = 3) +
    scale_color_manual(values = site_cols) +
    labs(x = "PC1 (54%)", y = "PC3 (12%)", color = "Site") +
    my_gg_theme.def_size +
    theme(
        legend.title = element_text(size = 12)
    )

p3 = ggplot() +
    geom_hline(yintercept = 0, lty = 2) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_segment(data = pca.spp, 
        aes(x = rep(0, nrow(pca.spp)), y = rep(0, nrow(pca.spp)), xend = MDS2, yend = MDS3),
        arrow = arrow(length = unit(0.1, "cm"))
    ) +
    geom_text(data = pca.spp, aes(x = MDS2*1.1, y = MDS3*1.1, label = var)) +
    geom_point(data = pca.sites, aes(x = MDS2, y = MDS3, color = state), size = 3) +
    scale_color_manual(values = site_cols) +
    labs(x = "PC2 (21%)", y = "PC3 (12%)", color = "Site") +
    my_gg_theme.def_size +
    theme(
        legend.title = element_text(size = 12)
    )

pdf("figures/GxE/climate_vars_pca.pdf", width = 7.5, height = 6)
p1
p2
p3
dev.off()

######################################
######################################
#pairwise cors
site_metadata.vars = site_metadata %>% 
        select(HDD4.mean_growing, HDD4.mean_nongrowing, tmin, tmax, MAT, ppt, freezeThaw.annual_mean, duration_infection)

cor(site_metadata.vars)
pairs(site_metadata.vars)

site_metadata.vars = site_metadata %>% 
        select(state, HDD4.mean_growing, HDD4.mean_nongrowing, tmin, tmax, MAT, ppt, freezeThaw.annual_mean, duration_infection)

colnames(site_metadata.vars)
col_labs = c(
    #"state",
    "GDD (growing)",  
    "GDD (nongrowing)", 
    "Tmin",
    "Tmax", 
    "MAT", 
    "MAP", 
    "freeze-thaw",
    "Infestation age"
)

p1 = ggpairs(
    site_metadata.vars, 
    columns = 2:ncol(site_metadata.vars),
    columnLabels = col_labs
) +
    my_gg_theme.def_size +
    theme(
        axis.text.x = element_text(angle = 55, hjust = 1)
        )
p1

pdf("figures/GxE/climate_vars_pairs_cor.pdf", width = 10, height = 10)
p1
dev.off()
