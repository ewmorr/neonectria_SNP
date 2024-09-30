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


#convert to nt difs per kb
dist.Nf = dist.Nf / (41018940 / 1000)
dist.Nd = dist.Nd / (38535154 / 1000)

#get order of samples
dist.order.Nf = as.matrix(dist.Nf ) %>% rownames
dist.order.Nd = as.matrix(dist.Nd ) %>% rownames 

length(dist.order.Nf)
length(dist.order.Nd)

#metadata
sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
sample_metadata.Nd = read.csv("data/sample_metadata/Nd_filtered.lat_lon_dur_inf.csv")

nrow(sample_metadata.Nf)
nrow(sample_metadata.Nd)
rownames(sample_metadata.Nf) = sample_metadata.Nf$Sequence_label
rownames(sample_metadata.Nd) = sample_metadata.Nd$Sequence_label

#lat lon tables (ordered)
dur.Nf = data.frame(
    samp = sample_metadata.Nf[dist.order.Nf, "Sequence_label"],
    collection_period = sample_metadata.Nf[dist.order.Nf, "collection_period"],
    lat = sample_metadata.Nf[dist.order.Nf, "lat"],
    lon = sample_metadata.Nf[dist.order.Nf, "lon"],
    dur_inf = sample_metadata.Nf[dist.order.Nf, "duration_infection"],
    stringsAsFactors = T
)
rownames(dur.Nf) = dur.Nf$samp

dur.Nd = data.frame(
    samp = sample_metadata.Nd[dist.order.Nd, "Sequence_label"],
    collection_period = sample_metadata.Nd[dist.order.Nd, "collection_period"],
    lat = sample_metadata.Nd[dist.order.Nd, "lat"],
    lon = sample_metadata.Nd[dist.order.Nd, "lon"],
    dur_inf = sample_metadata.Nd[dist.order.Nd, "duration_infection"],
    stringsAsFactors = T
)
rownames(dur.Nd) = dur.Nd$samp


#For now we only want to compare isolates collected in the modern era
# We may use the older isolates for a different test
dur.Nf$lat[dur.Nf$collection_period == "early"] = NA
dur.Nf$lon[dur.Nf$collection_period == "early"] = NA
dur.Nf$dur_inf[dur.Nf$collection_period == "early"] = NA

dur.Nd$lat[dur.Nd$collection_period == "early"] = NA
dur.Nd$lon[dur.Nd$collection_period == "early"] = NA
dur.Nd$dur_inf[dur.Nd$collection_period == "early"] = NA

#calculate geographic distance
Nf.Dgeo <- distm(x = dur.Nf[,c("lon", "lat")], fun = distVincentyEllipsoid)
Nd.Dgeo <- distm(x = dur.Nd[,c("lon", "lat")], fun = distVincentyEllipsoid)
Nf.Ddur <- dist(x = dur.Nf$dur_inf) %>% as.matrix
Nd.Ddur <- dist(x = dur.Nd$dur_inf) %>% as.matrix

rownames(Nf.Ddur) = dur.Nf$samp
colnames(Nf.Ddur) = dur.Nf$samp

rownames(Nd.Ddur) = dur.Nd$samp
colnames(Nd.Ddur) = dur.Nd$samp

#set duration difference to 0 based on zero distances to NA (within site comps)
sum(Nf.Dgeo == 0, na.rm = T)
Nf.Ddur[Nf.Dgeo == 0] = NA
sum(Nd.Dgeo == 0, na.rm = T)
Nd.Ddur[Nd.Dgeo == 0] = NA
Nf.Ddur = as.dist(Nf.Ddur)
Nd.Ddur = as.dist(Nd.Ddur)

# match NAs in gen dists
dist.Nf[is.na(Nf.Ddur)] = NA
dist.Nd[is.na(Nd.Ddur)] = NA

sum(is.na(dist.Nf))
sum(is.na(Nf.Ddur))
sum(is.na(dist.Nd))
sum(is.na(Nd.Ddur))

mantel(dist.Nf, Nf.Ddur, na.rm = T)
#Mantel statistic r: 0.2144 
#Significance: 0.001 
mantel(dist.Nd, Nd.Ddur, na.rm = T)
#Mantel statistic r: 0.01163 
#Significance: 0.444 


#rm VA from Nf and rerun
VA_samps = sample_metadata.Nf %>% filter(state == "VA") %>% pull(Sequence_label)
dist.Nf.no_VA = as.matrix(dist.Nf)
dist.Nf.no_VA[rownames(as.matrix(dist.Nf)) %in% VA_samps, ] = NA
dist.Nf.no_VA[,colnames(as.matrix(dist.Nf)) %in% VA_samps] = NA

Nf.Ddur.no_VA = as.matrix(Nf.Ddur)
Nf.Ddur.no_VA[is.na(dist.Nf.no_VA)] = NA

dist.Nf.no_VA = as.dist(dist.Nf.no_VA)
Nf.Ddur.no_VA = as.dist(Nf.Ddur.no_VA)

mantel(dist.Nf.no_VA, Nf.Ddur.no_VA, na.rm = T)
#Mantel statistic r: 0.1599 
#Significance: 0.001 


#########################
#long format for plotting
Nf.Ddur.long = reshape2::melt(Nf.Ddur %>% as.matrix)
Nf.Ddur.long %>% filter(Var1 == "NG2" & Var2 == "NG1") %>% nrow
Nd.Ddur.long = reshape2::melt(Nd.Ddur %>% as.matrix)
#Need to set self comps to NA (zero dists, i.e., same site, have been removed already..)
sum(is.na(Nf.Ddur.long$value))
Nf.Ddur.long[Nf.Ddur.long$Var1 == Nf.Ddur.long$Var2, "value"] = NA
sum(is.na(Nf.Ddur.long$value))
Nd.Ddur.long[Nd.Ddur.long$Var1 == Nd.Ddur.long$Var2, "value"] = NA

sum(is.na(Nf.Ddur.long$value))
sum(is.na(Nf.Dgen.long$value))
#set gen dists to NA where dur dist == NA
Nf.Dgen.long = reshape2::melt(dist.Nf %>% as.matrix)
Nd.Dgen.long = reshape2::melt(dist.Nd %>% as.matrix)

Nf.Dgen.long[is.na(Nf.Ddur.long$value), "value"] = NA
sum(is.na(Nf.Ddur.long$value))
sum(is.na(Nf.Dgen.long$value))
Nd.Dgen.long[is.na(Nd.Ddur.long$value), "value"] = NA

colnames(Nf.Ddur.long)[3] = "durDif"
colnames(Nf.Dgen.long)[3] = "SNPsPerKb"
colnames(Nd.Ddur.long)[3] = "durDif"
colnames(Nd.Dgen.long)[3] = "SNPsPerKb"

Nf.long = full_join(Nf.Ddur.long, Nf.Dgen.long)
nrow(Nf.long)
nrow(Nf.Ddur.long)
nrow(Nf.Dgen.long)
Nd.long = full_join(Nd.Ddur.long, Nd.Dgen.long)

######################
#plot
#

cor.test(Nf.long$durDif, Nf.long$SNPsPerKb, na.action = na.rm)

p1 = ggplot(Nf.long, aes(x = durDif, y = SNPsPerKb)) +
    geom_point(alpha = 0.08, shape = 1) +
    geom_smooth(method = "lm", linetype = 2, color = "black") +
    labs(x = "Difference in infestation duration (years)", y = "Hamming distance (SNPs per Kb)", title = "a") +
    annotate(
        geom = "text", 
        label = expression(paste("Mantel r = 0.21, ", italic("P"), " = 0.001")),
        x = 61,
        y = 2.5
    ) +
    scale_y_continuous(breaks = c(2,3,4,5)) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.10)
    )
p1

p2 = ggplot(Nd.long, aes(x = durDif, y = SNPsPerKb)) +
    geom_point(alpha = 0.25, shape = 1) +
    geom_smooth(method = "lm", linetype = 2, color = "black") +
    labs(x = "Difference in infestation duration (years)", y = "Hamming distance (SNPs per Kb)", title = "b") +
    scale_y_continuous(breaks = c(1,4,7,10)) +
    annotate(
        geom = "text", 
        label = expression(paste("Mantel r = 0.01, ", italic("P"), " = 0.44")),
        x = 69,
        y = 0.75
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.06),
        axis.title.y = element_blank()
    )
p2

pdf("figures/pop_gen/IBD/IB_durationInfection.pdf", width = 10, height = 3.5)
grid.arrange(p1,p2,ncol = 2)
dev.off()

