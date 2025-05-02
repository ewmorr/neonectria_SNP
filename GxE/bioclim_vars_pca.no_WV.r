#read env data from sample list and unique
library(vegan)
library(dplyr)
library(ggplot2)
library(GGally)
source("library/ggplot_theme.txt")
library(FactoMineR)
library(factoextra)

#Join pops to site data
site.info = read.csv("data/sample_metadata/site_info.csv")

site.bioclim = read.csv("data/sample_metadata/Nf.site_bioclim.csv")
site.bioclim = site.bioclim %>% filter(state != "WV")

source("data/sample_metadata/world_clim/bioclim_names.txt")
bioclim_var_names

# can pull in dur inf
site_bioclim_w_dur_inf = left_join(
    site.bioclim,
    site.info %>% select(state, lat, lon, duration_infection),
    by = c("state", "lat", "lon")
) %>% unique()

nrow(site_bioclim_w_dur_inf)

#scale the vars
site_metadata.scaled = apply(
    site.bioclim %>% 
        select(-c(state, lat, lon)), 
    2, 
    scale
) %>% as.data.frame

rownames(site_metadata.scaled) = site.bioclim$state
head(site_metadata.scaled)

#run princomp
climate.pca = capscale(site_metadata.scaled ~ 1, distance = "euclidean")

str(climate.pca)
plot(climate.pca$CA$eig/sum(climate.pca$CA$eig) )
climate.pca$CA$eig/sum(climate.pca$CA$eig)
#         MDS1         MDS2         MDS3         MDS4         MDS5         MDS6         MDS7         MDS8         MDS9 
# 5.614439e-01 2.080532e-01 1.265160e-01 5.199811e-02 3.689990e-02 6.884275e-03 3.958608e-03 2.278237e-03 7.051400e-04 
#        MDS10        MDS11        MDS12        MDS13        MDS14        MDS15        MDS16        MDS17        MDS18 
# 5.951477e-04 2.606305e-04 1.837469e-04 1.068166e-04 6.169814e-05 4.301795e-05 1.028420e-05 1.250954e-06 7.103047e-08 

sum((climate.pca$CA$eig/sum(climate.pca$CA$eig))[1:2])
#First two axes contain  0.7694971
sum((climate.pca$CA$eig/sum(climate.pca$CA$eig))[1:3])
#First three axes contain  0.8960131
sum((climate.pca$CA$eig/sum(climate.pca$CA$eig))[1:5])
#First five axes contain 0.9849111 of var

biplot(climate.pca)
biplot(climate.pca, choices = c(1,3))
biplot(climate.pca, choices = c(2,3))

pdf("figures/GxE/bioclim_cors/PCA_w_NC.pdf")
biplot(climate.pca, xlab = "PC1 (56% variance)", ylab = "PC2 (21% variance)")
dev.off()


#####################################
#####################################
# remove NC and re-examine cors
site.bioclim = read.csv("data/sample_metadata/Nf.site_bioclim.csv")
site.bioclim = site.bioclim %>% filter(state != "WV" & state != "NC")
#scale the vars
site_metadata.scaled = apply(
    site.bioclim %>% 
        select(-c(state, lat, lon)), 
    2, 
    scale
) %>% as.data.frame

rownames(site_metadata.scaled) = site.bioclim$state
head(site_metadata.scaled)

#run princomp
climate.pca = capscale(site_metadata.scaled ~ 1, distance = "euclidean")

str(climate.pca)
plot(climate.pca$CA$eig/sum(climate.pca$CA$eig) )
climate.pca$CA$eig/sum(climate.pca$CA$eig)
#         MDS1         MDS2         MDS3         MDS4         MDS5         MDS6         MDS7         MDS8         MDS9        MDS10        MDS11        MDS12        MDS13 
# 4.643976e-01 2.678374e-01 1.344070e-01 7.203952e-02 3.749587e-02 1.025207e-02 5.761864e-03 4.125099e-03 1.494391e-03 1.047478e-03 4.561089e-04 3.276941e-04 2.398909e-04 
#        MDS14        MDS15        MDS16        MDS17        MDS18 
# 7.704976e-05 2.618406e-05 1.341076e-05 1.243700e-06 9.441154e-08 

sum((climate.pca$CA$eig/sum(climate.pca$CA$eig))[1:2])
#First two axes contain  0.732235
sum((climate.pca$CA$eig/sum(climate.pca$CA$eig))[1:3])
#First three axes contain  0.866642
sum((climate.pca$CA$eig/sum(climate.pca$CA$eig))[1:5])
#First five axes contain 0.9761774 of var

biplot(climate.pca)
biplot(climate.pca, choices = c(1,3))
biplot(climate.pca, choices = c(2,3))

pdf("figures/GxE/bioclim_cors/PCA_no_NC.pdf")
biplot(climate.pca, xlab = "PC1 (46% variance)", ylab = "PC2 (27% variance)")
dev.off()

######################################
######################################
######################################
# Compare to FactoMineR
######################################
######################################
site.bioclim = read.csv("data/sample_metadata/Nf.site_bioclim.csv")
site.bioclim = site.bioclim %>% filter(state != "WV")
site_metadata.scaled = apply(
    site.bioclim %>% 
        select(-c(state, lat, lon)), 
    2, 
    scale
) %>% as.data.frame

site_metadata.scaled.names = site_metadata.scaled
colnames(site_metadata.scaled.names) = bioclim_var_names

climate.factoPCA = FactoMineR::PCA(
    site_metadata.scaled, 
    scale.unit = F, #already scaled
    ncp = 5#the def is five, which is OK for this dataset, but may want more
)

climate.factoPCA = FactoMineR::PCA(
    site_metadata.scaled.names, 
    scale.unit = F, #already scaled
    ncp = 5#the def is five, which is OK for this dataset, but may want more
)

plot(climate.factoPCA, choix  = "var", axes = c(1,2))
plot(climate.factoPCA, choix  = "var", axes = c(1,3))

str(climate.factoPCA)

climate.factoPCA$var$contrib

# Contributions of variables to PC1
fviz_contrib(climate.factoPCA, choice = "var", axes = 1, top = 20)
# temp seasonality and precip; temp of driest (+), annual precip (+), temp seasonality (-), precip wettest (+), 
fviz_contrib(climate.factoPCA, choice = "ind", axes = 1, top = 24)

# Contributions of variables to PC2
fviz_contrib(climate.factoPCA, choice = "var", axes = 2, top = 20)   
# temp; temp warmest quarter and month, MAT
fviz_contrib(climate.factoPCA, choice = "ind", axes = 2, top = 24)

# Contributions of variables to PC3
fviz_contrib(climate.factoPCA, choice = "var", axes = 3, top = 20)    
# temp wettest and precip seasonality
fviz_contrib(climate.factoPCA, choice = "ind", axes = 3, top = 24)

# Contributions of variables to PC4
fviz_contrib(climate.factoPCA, choice = "var", axes = 4, top = 20)    
# this is very clearly a temp variatino axis; diurnal range is top by a lot
# then isothermality (diurnal range/annual range)
fviz_contrib(climate.factoPCA, choice = "ind", axes = 4, top = 24)

# Contributions of variables to PC5
fviz_contrib(climate.factoPCA, choice = "var", axes = 5, top = 20)    
# temp of wettest by far strongest
# then temp driest
fviz_contrib(climate.factoPCA, choice = "ind", axes = 5, top = 24)



######################################
######################################
#pairwise cors
site_metadata.scaled
cor(site_metadata.scaled)
pairs(site_metadata.scaled)

p1 = ggpairs(
    site_metadata.scaled, 
    columns = 1:ncol(site_metadata.scaled)#,
#    columnLabels = bioclim_var_names
) +
    my_gg_theme.def_size +
    theme(
        axis.text.x = element_text(angle = 55, hjust = 1)
        )
p1

pdf("figures/GxE/Nf_sites.bioclim.pairs_cor.no_WV.pdf", width = 16, height = 16)
p1
dev.off()


#################################
# running new PCAs after examine cors
# 
# 
# 

#we take the same sets of factors from below where we excluded NC. This is two
# sets of highly correlated vars
# all temps except for max of wettest and diurnal range
# all precips together plus max temp of wettest

# temps pca
climate.factoPCA.temps = FactoMineR::PCA(
    site_metadata.scaled.names[,c(1,3:7,9:11)], 
    scale.unit = F, #already scaled
    ncp = 5#the def is five, which is OK for this dataset, but may want more
)

# Contributions of variables to PC1
fviz_contrib(climate.factoPCA.temps, choice = "var", axes = 1, top = 20)
# cold temps and the amount of seasonal/annual variability
fviz_contrib(climate.factoPCA.temps, choice = "ind", axes = 1, top = 24)

# Contributions of variables to PC2
fviz_contrib(climate.factoPCA.temps, choice = "var", axes = 2, top = 20)   
# warm temps
fviz_contrib(climate.factoPCA.temps, choice = "ind", axes = 2, top = 24)
# PC2 is just max temp warmest and mean temp warmest


# precip pca
climate.factoPCA.precip = FactoMineR::PCA(
    site_metadata.scaled.names[,c(8,12:19)], 
    scale.unit = F, #already scaled
    ncp = 5#the def is five, which is OK for this dataset, but may want more
)

# Contributions of variables to PC1
fviz_contrib(climate.factoPCA.precip, choice = "var", axes = 1, top = 20)
# MAP and precip of driest,coldest, wettest
fviz_contrib(climate.factoPCA, choice = "ind", axes = 1, top = 24)

# Contributions of variables to PC2
fviz_contrib(climate.factoPCA.precip, choice = "var", axes = 2, top = 20)   
# seasonality
fviz_contrib(climate.factoPCA, choice = "ind", axes = 2, top = 24)
#temp wettest, precip seasonality, precip warmest
# seasonality axis

PC.df = data.frame(
    temp.PC1 = climate.factoPCA.temps$ind$coord[,1],
    #temp.PC2 = climate.factoPCA.temps$ind$coord[,2],
    precip.PC1 = climate.factoPCA.precip$ind$coord[,1],
    precip.PC2 = climate.factoPCA.precip$ind$coord[,2],
    bio2 = site_metadata.scaled.names[,2]
)

# look at cors between new vars
cor(PC.df)
pairs(PC.df)

p2 = ggpairs(
    PC.df, 
    columns = 1:ncol(PC.df)#,
#    columnLabels = bioclim_var_names
) +
    my_gg_theme.def_size +
    theme(
        axis.text.x = element_text(angle = 55, hjust = 1)
    )
p2

pdf("figures/GxE/Nf_sites.bioclim_PCA_temp_precip_separate.pairs_cor.no_WV.pdf", width = 16, height = 16)
p2
dev.off()


# we still end up with some fairly strong multicollinearity. let's stick with the full PCA approach
#  and we interpret the best we can

#PC.df = data.frame(
#    bioclim.PC1 = climate.factoPCA$ind$coord[,1],
#    bioclim.PC2 = climate.factoPCA$ind$coord[,2],
#    bioclim.PC3 = climate.factoPCA$ind$coord[,3],
#    bioclim.PC4 = climate.factoPCA$ind$coord[,4],
#    bioclim.PC5 = climate.factoPCA$ind$coord[,5]
#)

PC.df$state = site.bioclim$state

write.csv(PC.df, "data/sample_metadata/Nf.sites_bioclim_PC.no_WV.csv", row.names = F, quote = F)

climate.factoPCA$eig
#          eigenvalue percentage of variance cumulative percentage of variance
#comp 1  1.020363e+01           5.614439e+01                          56.14439
#comp 2  3.781140e+00           2.080532e+01                          76.94971
#comp 3  2.299291e+00           1.265160e+01                          89.60131
#comp 4  9.450092e-01           5.199811e+00                          94.80112
#comp 5  6.706156e-01           3.689990e+00                          98.49111
#comp 6  1.251142e-01           6.884275e-01                          99.17954
#comp 7  7.194340e-02           3.958608e-01                          99.57540
#comp 8  4.140448e-02           2.278237e-01                          99.80322
#comp 9  1.281515e-02           7.051400e-02                          99.87373
#comp 10 1.081616e-02           5.951477e-02                          99.93325
#comp 11 4.736677e-03           2.606305e-02                          99.95931
#comp 12 3.339400e-03           1.837469e-02                          99.97769
#comp 13 1.941276e-03           1.068166e-02                          99.98837
#comp 14 1.121297e-03           6.169814e-03                          99.99454
#comp 15 7.818044e-04           4.301795e-03                          99.99884
#comp 16 1.869041e-04           1.028420e-03                          99.99987
#comp 17 2.273473e-05           1.250954e-04                          99.99999
#comp 18 1.290902e-06           7.103047e-06                         100.00000
#comp 19 1.681441e-29           9.251950e-29                         100.00000



######################################
######################################
######################################
# Re-running the whole thing without NC
######################################
######################################
site.bioclim = read.csv("data/sample_metadata/Nf.site_bioclim.csv")
site.bioclim = site.bioclim %>% filter(state != "WV" & state != "NC")
site_metadata.scaled = apply(
    site.bioclim %>% 
        select(-c(state, lat, lon)), 
    2, 
    scale
) %>% as.data.frame

site_metadata.scaled.names = site_metadata.scaled
colnames(site_metadata.scaled.names) = bioclim_var_names

climate.factoPCA = FactoMineR::PCA(
    site_metadata.scaled, 
    scale.unit = F, #already scaled
    ncp = 5#the def is five, which is OK for this dataset, but may want more
)

climate.factoPCA = FactoMineR::PCA(
    site_metadata.scaled.names, 
    scale.unit = F, #already scaled
    ncp = 5#the def is five, which is OK for this dataset, but may want more
)
# we now get two groups of opposite loadings
# mean temp and other temp variables *except* seasonality load + 1 and - 2
# temp seasonality load - 1 and + 2
# mean precip and all except seasonality load + 1 and + 2
# precip seasonality load - 1 and - 2

plot(climate.factoPCA, choix  = "var", axes = c(1,2))
plot(climate.factoPCA, choix  = "var", axes = c(1,3))

str(climate.factoPCA)

climate.factoPCA$var$contrib

# Contributions of variables to PC1
fviz_contrib(climate.factoPCA, choice = "var", axes = 1, top = 20)
# temp seasonality and precip; temp of driest (+), annual precip (+), temp seasonality (-), precip wettest (+), 
fviz_contrib(climate.factoPCA, choice = "ind", axes = 1, top = 24)

# Contributions of variables to PC2
fviz_contrib(climate.factoPCA, choice = "var", axes = 2, top = 20)   
# temp; temp warmest quarter and month, MAT
fviz_contrib(climate.factoPCA, choice = "ind", axes = 2, top = 24)

# Contributions of variables to PC3
fviz_contrib(climate.factoPCA, choice = "var", axes = 3, top = 20)    
# temp wettest and precip seasonality
fviz_contrib(climate.factoPCA, choice = "ind", axes = 3, top = 24)

# Contributions of variables to PC4
fviz_contrib(climate.factoPCA, choice = "var", axes = 4, top = 20)    
# this is very clearly a temp variatino axis; diurnal range is top by a lot
# then isothermality (diurnal range/annual range)
fviz_contrib(climate.factoPCA, choice = "ind", axes = 4, top = 24)

# Contributions of variables to PC5
fviz_contrib(climate.factoPCA, choice = "var", axes = 5, top = 20)    
# temp of wettest by far strongest
# then temp driest
fviz_contrib(climate.factoPCA, choice = "ind", axes = 5, top = 24)



######################################
######################################
#pairwise cors
site_metadata.scaled
cor(site_metadata.scaled)
pairs(site_metadata.scaled)

p1 = ggpairs(
    site_metadata.scaled, 
    columns = 1:ncol(site_metadata.scaled)#,
#    columnLabels = bioclim_var_names
) +
    my_gg_theme.def_size +
    theme(
        axis.text.x = element_text(angle = 55, hjust = 1)
        )
p1

pdf("figures/GxE/bioclim_cors/Nf_sites.bioclim.pairs_cor.no_WV_no_NC.pdf", width = 16, height = 16)
p1
dev.off()


#################################
# running new PCAs after examine cors
# 
# 
# 

#keeping mean temp wettest out becasue it is highly ocrreltated wiuth precip seasonality
site_metadata.scaled.names

# temps pca
climate.factoPCA.temps = FactoMineR::PCA(
    site_metadata.scaled.names[,c(1,3:7,9:11)], 
    scale.unit = F, #already scaled
    ncp = 5#the def is five, which is OK for this dataset, but may want more
)
# positive higher temp, negative more seasonality/variation
# #81.39% variance.
# only need axis 1

# 
# Contributions of variables to PC1
fviz_contrib(climate.factoPCA.temps, choice = "var", axes = 1, top = 20)
# cold temps and the amount of seasonal/annual variability
fviz_contrib(climate.factoPCA.temps, choice = "ind", axes = 1, top = 24)

# Contributions of variables to PC2
fviz_contrib(climate.factoPCA.temps, choice = "var", axes = 2, top = 20)   
# warm temps
fviz_contrib(climate.factoPCA.temps, choice = "ind", axes = 2, top = 24)


# precip pca
climate.factoPCA.precip = FactoMineR::PCA(
    site_metadata.scaled.names[,c(8,12:19)], 
    scale.unit = F, #already scaled
    ncp = 5#the def is five, which is OK for this dataset, but may want more
)
# axis 1 positive more precip, negative more seasonality
# 69.01% variance
# axis 2 [positive warm and wet (temp wettest and precip warmest)
# 20.2% variance

# Contributions of variables to PC1
fviz_contrib(climate.factoPCA.precip, choice = "var", axes = 1, top = 20)
# MAP and precip of driest,coldest, wettest
fviz_contrib(climate.factoPCA, choice = "ind", axes = 1, top = 24)

# Contributions of variables to PC2
fviz_contrib(climate.factoPCA.precip, choice = "var", axes = 2, top = 20)   
# seasonality
fviz_contrib(climate.factoPCA, choice = "ind", axes = 2, top = 24)


PC.df = data.frame(
    temp.PC1 = climate.factoPCA.temps$ind$coord[,1],
    #temp.PC2 = climate.factoPCA.temps$ind$coord[,2],
    precip.PC1 = climate.factoPCA.precip$ind$coord[,1],
    precip.PC2 = climate.factoPCA.precip$ind$coord[,2],
    bio2 = site_metadata.scaled.names[,2]
)

# look at cors between new vars
cor(PC.df)
pairs(PC.df)

p2 = ggpairs(
    PC.df, 
    columns = 1:ncol(PC.df)#,
#    columnLabels = bioclim_var_names
) +
    my_gg_theme.def_size +
    theme(
        axis.text.x = element_text(angle = 55, hjust = 1)
    )
p2

pdf("figures/GxE/Nf_sites.bioclim_PCA.pairs_cor.no_WV_no_NC.pdf", width = 16, height = 16)
p1
dev.off()



write.csv(PC.df, "data/sample_metadata/Nf.sites_bioclim_PC.no_WV_no_NC.csv", row.names = F, quote = F)


climate.factoPCA.temps$eig
#         eigenvalue percentage of variance cumulative percentage of variance
#comp 1 6.992405e+00           8.139308e+01                          81.39308
#comp 2 1.045852e+00           1.217394e+01                          93.56702
#comp 3 2.774816e-01           3.229944e+00                          96.79696
#comp 4 2.488145e-01           2.896254e+00                          99.69321
#comp 5 2.160732e-02           2.515137e-01                          99.94473
#comp 6 4.119886e-03           4.795635e-02                          99.99268
#comp 7 5.484718e-04           6.384328e-03                          99.99907
#comp 8 8.007244e-05           9.320601e-04                         100.00000
#comp 9 7.757803e-29           9.030247e-28                         100.00000

climate.factoPCA.precip$eig
#        eigenvalue percentage of variance cumulative percentage of variance
#comp 1 5.928394173            69.00776286                          69.00776
#comp 2 1.734004973            20.18418487                          89.19195
#comp 3 0.689243363             8.02293861                          97.21489
#comp 4 0.158108989             1.84042210                          99.05531
#comp 5 0.033810325             0.39355933                          99.44887
#comp 6 0.029051559             0.33816629                          99.78703
#comp 7 0.010061051             0.11711277                          99.90415
#comp 8 0.005102617             0.05939555                          99.96354
#comp 9 0.003132041             0.03645762                         100.00000


