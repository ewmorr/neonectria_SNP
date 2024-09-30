library(dplyr)
library(ggplot2)
library(gridExtra)
library(vegan)
library(geosphere)
source("library/ggplot_theme.txt")

dist.Nf = read.table("data/Nf/final_tables/rm_dups/FINAL_invariant.IBD_analyses.dist", header = F) 
dist.ID.Nf = read.table("data/Nf/final_tables/rm_dups/FINAL_invariant.IBD_analyses.dist.id", header = F)
rownames(dist.Nf) = dist.ID.Nf[,1]
colnames(dist.Nf) = dist.ID.Nf[,2]
dist.Nf = dist.Nf %>% as.dist()
dist.Nf

dist.Nd = read.table("data/Nd/final_tables/rm_dups/FINAL_invariant.IBD_analyses.dist", header = F) 
dist.ID.Nd = read.table("data/Nd/final_tables/rm_dups/FINAL_invariant.IBD_analyses.dist.id", header = F)
rownames(dist.Nd) = dist.ID.Nd[,1]
colnames(dist.Nd) = dist.ID.Nd[,2]
dist.Nd = dist.Nd %>% as.dist()
head(dist.Nd)

dist.Nc = read.table("data/Nc/final_tables/FINAL_invariant.IBD_analyses.dist", header = F) 
dist.ID.Nc = read.table("data/Nc/final_tables/FINAL_invariant.IBD_analyses.dist.id", header = F)
rownames(dist.Nc) = dist.ID.Nc[,1]
colnames(dist.Nc) = dist.ID.Nc[,2]
dist.Nc = dist.Nc %>% as.dist()
head(dist.Nc)

#convert to nt difs per kb
dist.Nf = dist.Nf / (41018940 / 1000)
dist.Nd = dist.Nd / (38535154 / 1000)
dist.Nc = dist.Nc / (40630626 / 1000)

#get order of samples
dist.order.Nf = as.matrix(dist.Nf ) %>% rownames
dist.order.Nd = as.matrix(dist.Nd ) %>% rownames 
dist.order.Nc = as.matrix(dist.Nc ) %>% rownames

length(dist.order.Nf)
length(dist.order.Nd)
length(dist.order.Nc)

#metadata
sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
sample_metadata.Nd = read.csv("data/sample_metadata/Nd_filtered.lat_lon_dur_inf.csv")
sample_metadata.Nc = read.csv("data/sample_metadata/Nc_canton_loc_date.lat_lon.csv")

nrow(sample_metadata.Nf)
nrow(sample_metadata.Nd)
nrow(sample_metadata.Nc)
rownames(sample_metadata.Nf) = sample_metadata.Nf$Sequence_label
rownames(sample_metadata.Nd) = sample_metadata.Nd$Sequence_label
rownames(sample_metadata.Nc) = sample_metadata.Nc$Sequence_label

#lat lon tables (ordered)
coords.Nf = data.frame(
    samp = sample_metadata.Nf[dist.order.Nf, "Sequence_label"],
    collection_period = sample_metadata.Nf[dist.order.Nf, "collection_period"],
    lat = sample_metadata.Nf[dist.order.Nf, "lat"],
    lon = sample_metadata.Nf[dist.order.Nf, "lon"],
    stringsAsFactors = T
)
rownames(coords.Nf) = coords.Nf$samp

coords.Nd = data.frame(
    samp = sample_metadata.Nd[dist.order.Nd, "Sequence_label"],
    collection_period = sample_metadata.Nd[dist.order.Nd, "collection_period"],
    lat = sample_metadata.Nd[dist.order.Nd, "lat"],
    lon = sample_metadata.Nd[dist.order.Nd, "lon"],
    stringsAsFactors = T
)
rownames(coords.Nd) = coords.Nd$samp

coords.Nc = data.frame(
    samp = sample_metadata.Nc[dist.order.Nc, "Sequence_label"],
    collection_period = rep("modern", nrow(sample_metadata.Nc)),
    lat = sample_metadata.Nc[dist.order.Nc, "lat"],
    lon = sample_metadata.Nc[dist.order.Nc, "lon"],
    stringsAsFactors = T
)
rownames(coords.Nc) = coords.Nc$samp

#For now we only want to compare isolates collected in the modern era
# We may use the older isolates for a different test
coords.Nf$lat[coords.Nf$collection_period == "early"] = NA
coords.Nf$lon[coords.Nf$collection_period == "early"] = NA

coords.Nd$lat[coords.Nd$collection_period == "early"] = NA
coords.Nd$lon[coords.Nd$collection_period == "early"] = NA

#calculate geographic distance
Nf.Dgeo <- distm(x = coords.Nf[,c("lon", "lat")], fun = distVincentyEllipsoid)
Nd.Dgeo <- distm(x = coords.Nd[,c("lon", "lat")], fun = distVincentyEllipsoid)
Nc.Dgeo <- distm(x = coords.Nc[,c("lon", "lat")], fun = distVincentyEllipsoid)

rownames(Nf.Dgeo) = coords.Nf$samp
colnames(Nf.Dgeo) = coords.Nf$samp
Nf.Dgeo = as.dist(Nf.Dgeo/1000)

rownames(Nd.Dgeo) = coords.Nd$samp
colnames(Nd.Dgeo) = coords.Nd$samp
Nd.Dgeo = as.dist(Nd.Dgeo/1000)

rownames(Nc.Dgeo) = coords.Nc$samp
colnames(Nc.Dgeo) = coords.Nc$samp
Nc.Dgeo = as.dist(Nc.Dgeo/1000)


#set zero dists to NA (within site comps)
sum(Nf.Dgeo == 0, na.rm = T)
Nf.Dgeo[Nf.Dgeo == 0] = NA
sum(Nd.Dgeo == 0, na.rm = T)
Nd.Dgeo[Nd.Dgeo == 0] = NA
sum(Nc.Dgeo == 0, na.rm = T)
Nc.Dgeo[Nc.Dgeo == 0] = NA

# match NAs in gen dists
dist.Nf[is.na(Nf.Dgeo)] = NA
dist.Nd[is.na(Nd.Dgeo)] = NA
dist.Nc[is.na(Nc.Dgeo)] = NA


mantel(dist.Nf, Nf.Dgeo, na.rm = T)
#Mantel statistic r: 0.2465 
#Significance: 0.001 
mantel(dist.Nd, Nd.Dgeo, na.rm = T)
#Mantel statistic r: -0.08981 
#Significance: 0.782 
mantel(dist.Nc, Nc.Dgeo, na.rm = T)
#Mantel statistic r: 0.7574 
#Significance: 0.05 # restricted set of permutation

mantel(dist.Nf, log(Nf.Dgeo), na.rm = T)
#Mantel statistic r: 0.2307 
#Significance: 0.001 
mantel(dist.Nd, log(Nd.Dgeo), na.rm = T)
#Mantel statistic r: -0.05315 
#Significance: 0.675 
mantel(dist.Nc, log(Nc.Dgeo), na.rm = T)
#Mantel statistic r: 0.7849 
#Significance: 0.05 

#rm VA from Nf and rerun
VA_samps = sample_metadata.Nf %>% filter(state == "VA") %>% pull(Sequence_label)
dist.Nf.no_VA = as.matrix(dist.Nf)
dist.Nf.no_VA[rownames(as.matrix(dist.Nf)) %in% VA_samps, ] = NA
dist.Nf.no_VA[,colnames(as.matrix(dist.Nf)) %in% VA_samps] = NA

Nf.Dgeo.no_VA = as.matrix(Nf.Dgeo)
Nf.Dgeo.no_VA[is.na(dist.Nf.no_VA)] = NA

dist.Nf.no_VA = as.dist(dist.Nf.no_VA)
Nf.Dgeo.no_VA = as.dist(Nf.Dgeo.no_VA)

mantel(dist.Nf.no_VA, Nf.Dgeo.no_VA, na.rm = T)

#########################
#long format for plotting
Nf.Dgeo.long = reshape2::melt(Nf.Dgeo %>% as.matrix)
Nf.Dgeo.long %>% filter(Var1 == "NG2" & Var2 == "NG1") %>% nrow
Nd.Dgeo.long = reshape2::melt(Nd.Dgeo %>% as.matrix)
Nc.Dgeo.long = reshape2::melt(Nc.Dgeo %>% as.matrix)
#Need to set self comps to NA (zero dists, i.e., same site, have been removed already..)
Nf.Dgeo.long[Nf.Dgeo.long$Var1 == Nf.Dgeo.long$Var2, "value"] = NA
Nd.Dgeo.long[Nd.Dgeo.long$Var1 == Nd.Dgeo.long$Var2, "value"] = NA
Nc.Dgeo.long[Nc.Dgeo.long$Var1 == Nc.Dgeo.long$Var2, "value"] = NA

#set gen dists to NA where geo dist == NA
Nf.Dgen.long = reshape2::melt(dist.Nf %>% as.matrix)
Nd.Dgen.long = reshape2::melt(dist.Nd %>% as.matrix)
Nc.Dgen.long = reshape2::melt(dist.Nc %>% as.matrix)

Nf.Dgen.long[is.na(Nf.Dgeo.long$value), "value"] = NA
Nd.Dgen.long[is.na(Nd.Dgeo.long$value), "value"] = NA
Nc.Dgen.long[is.na(Nc.Dgeo.long$value), "value"] = NA

colnames(Nf.Dgeo.long)[3] = "km"
colnames(Nf.Dgen.long)[3] = "SNPsPerKb"
colnames(Nd.Dgeo.long)[3] = "km"
colnames(Nd.Dgen.long)[3] = "SNPsPerKb"
colnames(Nc.Dgeo.long)[3] = "km"
colnames(Nc.Dgen.long)[3] = "SNPsPerKb"

Nf.long = full_join(Nf.Dgeo.long, Nf.Dgen.long)
nrow(Nf.long)
nrow(Nf.Dgeo.long)
nrow(Nf.Dgen.long)
Nd.long = full_join(Nd.Dgeo.long, Nd.Dgen.long)
Nc.long = full_join(Nc.Dgeo.long, Nc.Dgen.long)

######################
#plot
#

cor.test(Nf.long$km, Nf.long$SNPsPerKb, na.action = na.rm)

p1 = ggplot(Nf.long, aes(x = km, y = SNPsPerKb)) +
    geom_point(alpha = 0.08, shape = 1) +
    geom_smooth(method = "lm", linetype = 2, color = "black") +
    scale_y_continuous(breaks = c(3,4,5)) +
    labs(x = "Geographic distance (km)", y = "Hamming distance (SNPs per Kb)", title = "a") +
    annotate(
        geom = "text", 
        label = expression(paste("Mantel r = 0.25, ", italic("P"), " = 0.001")),
        x = 1475,
        y = 2.5
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.10)
    )

p2 = ggplot(Nd.long, aes(x = km, y = SNPsPerKb)) +
    geom_point(alpha = 0.25, shape = 1) +
    geom_smooth(method = "lm", linetype = 2, color = "black") +
    labs(x = "Geographic distance (km)", y = "Hamming distance (SNPs per Kb)", title = "b") +
    scale_y_continuous(breaks = c(1,4,7,10)) +
    annotate(
        geom = "text", 
        label = expression(paste("Mantel r = -0.09, ", italic("P"), " = 0.79")),
        x = 1600,
        y = 0.75
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.06),
        axis.title.y = element_blank()
    )

p3 = ggplot(Nc.long, aes(x = km, y = SNPsPerKb)) +
    geom_point(alpha = 0.25, shape = 1) +
    geom_smooth(method = "lm", linetype = 2, color = "black") +
    labs(x = "Geographic distance (km)", y = "Hamming distance (SNPs per Kb)", title = "c") +
    annotate(
        geom = "text", 
        label = expression(paste("Mantel r = 0.75, ", italic("P"), " = 0.05")),
        x = 170,
        y = 12.5
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.12, margin = margin(b = -10))
    )

pdf("figures/pop_gen/IBD/IBD.pdf", width = 10, height = 3.5)
grid.arrange(p1,p2,ncol = 2)
dev.off()

pdf("figures/pop_gen/IBD/IBD_Nc.pdf", width = 6.5, height = 3.5)
p3
dev.off()
