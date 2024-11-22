library(dplyr)
library(ggplot2)
library(gridExtra)
library(vegan)
library(geosphere)
library(rlang)
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
rownames(Nf.Dgeo) = dur.Nf$samp
colnames(Nf.Dgeo) = dur.Nf$samp
Nf.Dgeo = as.dist(Nf.Dgeo/1000)

rownames(Nd.Ddur) = dur.Nd$samp
colnames(Nd.Ddur) = dur.Nd$samp
rownames(Nd.Dgeo) = dur.Nd$samp
colnames(Nd.Dgeo) = dur.Nd$samp
Nd.Dgeo = as.dist(Nd.Dgeo/1000)

#set zero dists to NA (within site comps)
sum(Nf.Dgeo == 0, na.rm = T)
#341
Nf.Dgeo[Nf.Dgeo == 0] = NA

sum(Nd.Dgeo == 0, na.rm = T)
#20
Nd.Dgeo[Nd.Dgeo == 0] = NA

Nf.Ddur = as.dist(Nf.Ddur)
Nd.Ddur = as.dist(Nd.Ddur)


# match NAs in gen dists and durs
dist.Nf[is.na(Nf.Dgeo)] = NA
Nf.Ddur[is.na(Nf.Dgeo)] = NA
dist.Nd[is.na(Nd.Dgeo)] = NA
Nd.Ddur[is.na(Nd.Dgeo)] = NA

sum(is.na(dist.Nf))
sum(is.na(Nf.Ddur))
sum(is.na(Nf.Dgeo))
sum(is.na(dist.Nd))
sum(is.na(Nd.Ddur))
sum(is.na(Nd.Dgeo))

####################
#Mantel geo dist
mantel(dist.Nf, Nf.Dgeo, na.rm = T)
#Mantel statistic r: 0.2467
#Significance: 0.001
mantel(dist.Nd, Nd.Dgeo, na.rm = T)
#Mantel statistic r: -0.09242
#Significance: 0.788

mantel(dist.Nf, log(Nf.Dgeo), na.rm = T)
#Mantel statistic r: 0.2304
#Significance: 0.001
mantel(dist.Nd, log(Nd.Dgeo), na.rm = T)
#Mantel statistic r: -0.04788
#Significance: 0.628

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
#Mantel statistic r: 0.2382 
#Significance: 0.002 

####################
#Mantel duration dist
mantel(dist.Nf, Nf.Ddur, na.rm = T)
#Mantel statistic r: 0.2144 
#Significance: 0.001 
mantel(dist.Nd, Nd.Ddur, na.rm = T)
#Mantel statistic r: 0.01163 
#Significance: 0.454 

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

#geo
Nf.Dgeo.mat = Nf.Dgeo %>% as.matrix
diag(Nf.Dgeo.mat) = NA
Nf.Dgeo.mat[upper.tri(Nf.Dgeo.mat)] = NA
Nf.Dgeo.long = reshape2::melt(Nf.Dgeo.mat %>% as.matrix)
Nf.Dgeo.long %>% filter(Var1 == "NG2" & Var2 == "NG1") %>% nrow
Nd.Dgeo.mat = Nd.Dgeo %>% as.matrix
diag(Nd.Dgeo.mat) = NA
Nd.Dgeo.mat[upper.tri(Nd.Dgeo.mat)] = NA
Nd.Dgeo.long = reshape2::melt(Nd.Dgeo %>% as.matrix)
#Need to set self comps to NA (zero dists, i.e., same site, have been removed already..)
##actually this is now handled above by diag() = NA
#Nf.Dgeo.long[Nf.Dgeo.long$Var1 == Nf.Dgeo.long$Var2, "value"] = NA
#Nd.Dgeo.long[Nd.Dgeo.long$Var1 == Nd.Dgeo.long$Var2, "value"] = NA

#set gen dists to NA where geo dist == NA
Nf.Dgen.mat = dist.Nf %>% as.matrix
diag(Nf.Dgen.mat) = NA
Nf.Dgen.mat[upper.tri(Nf.Dgen.mat)] = NA
Nf.Dgen.long = reshape2::melt(dist.Nf %>% as.matrix)
Nd.Dgen.mat = dist.Nd %>% as.matrix
diag(Nd.Dgen.mat) = NA
Nd.Dgen.mat[upper.tri(Nd.Dgen.mat)] = NA
Nd.Dgen.long = reshape2::melt(dist.Nd %>% as.matrix)

Nf.Dgen.long[is.na(Nf.Dgeo.long$value), "value"] = NA
Nd.Dgen.long[is.na(Nd.Dgeo.long$value), "value"] = NA

#duration
Nf.Ddur.mat = Nf.Ddur %>% as.matrix
diag(Nf.Ddur.mat) = NA
Nf.Ddur.mat[upper.tri(Nf.Ddur.mat)] = NA
Nf.Ddur.long = reshape2::melt(Nf.Ddur.mat)
Nf.Ddur.long %>% filter(Var1 == "NG121" & Var2 == "NG4") 
Nf.Ddur.long %>% filter(Var1 == "NG4" & Var2 == "NG121")
Nf.Dgeo.long %>% filter(Var1 == "NG121" & Var2 == "NG4") 
Nf.Dgeo.long %>% filter(Var1 == "NG4" & Var2 == "NG121")

Nd.Ddur.mat = Nd.Ddur %>% as.matrix
diag(Nd.Ddur.mat) = NA
Nd.Ddur.mat[upper.tri(Nd.Ddur.mat)] = NA
Nd.Ddur.long = reshape2::melt(Nd.Ddur.mat %>% as.matrix)
#Need to set self comps to NA (zero dists, i.e., same site, have been removed already..)
sum(is.na(Nf.Ddur.long$value))
Nf.Ddur.long[Nf.Ddur.long$Var1 == Nf.Ddur.long$Var2, "value"] = NA
sum(is.na(Nf.Ddur.long$value))
Nd.Ddur.long[Nd.Ddur.long$Var1 == Nd.Ddur.long$Var2, "value"] = NA

colnames(Nf.Dgeo.long)[3] = "km"
colnames(Nf.Ddur.long)[3] = "durDif"
colnames(Nf.Dgen.long)[3] = "SNPsPerKb"
colnames(Nd.Dgeo.long)[3] = "km"
colnames(Nd.Ddur.long)[3] = "durDif"
colnames(Nd.Dgen.long)[3] = "SNPsPerKb"

Nf.long = full_join(Nf.Dgeo.long, Nf.Ddur.long) %>%
    full_join(., Nf.Dgen.long)
nrow(Nf.long)
nrow(Nf.Dgeo.long)
nrow(Nf.Dgen.long)
Nd.long = full_join(Nd.Dgeo.long, Nd.Ddur.long) %>%
    full_join(., Nd.Dgen.long)

sum(is.na(Nf.Dgeo.long$km))
sum(is.na(Nf.Ddur.long$durDif))
sum(is.na(Nf.Dgen.long$SNPsPerKb))


######################
#plot
#

Nf.geo.mantel = mantel(dist.Nf, Nf.Dgeo, na.rm = T)
Nd.geo.mantel = mantel(dist.Nd, Nd.Dgeo, na.rm = T)
Nf.dur.mantel = mantel(dist.Nf, Nf.Ddur, na.rm = T)
Nd.dur.mantel = mantel(dist.Nd, Nd.Ddur, na.rm = T)


cor.test(Nf.long$km, Nf.long$SNPsPerKb, na.action = na.rm)
cor.test(Nf.long$durDif, Nf.long$SNPsPerKb, na.action = na.rm)

cor.test(Nd.long$km, Nd.long$SNPsPerKb, na.action = na.rm)
cor.test(Nd.long$durDif, Nd.long$SNPsPerKb, na.action = na.rm)

#these match the Mantel test statistic

Nf.geo.mantel.stat = round(Nf.geo.mantel$statistic,2)
Nf.geo.mantel.sig = round(Nf.geo.mantel$signif,3)

p1 = ggplot(Nf.long, aes(x = km, y = SNPsPerKb)) +
    geom_point(alpha = 0.08, shape = 1) +
    geom_smooth(method = "lm", linetype = 1, color = "white", linewidth = 2) +
    geom_smooth(method = "lm", linetype = 1, color = "black") +
    scale_y_continuous(breaks = c(3,4,5)) +
    labs(x = "Geographic distance (km)", y = "Hamming distance (SNPs per Kb)", title = "a") +
    annotate(
        geom = "text", 
        label = expr(paste("Mantel r = ", !!Nf.geo.mantel.stat, ", ", italic("P"), " = ", !!Nf.geo.mantel.sig)),
        x = 1475,
        y = 2.5
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.10)
    )
p1

Nd.geo.mantel.stat = round(Nd.geo.mantel$statistic,2)
Nd.geo.mantel.sig = round(Nd.geo.mantel$signif,2)

p2 = ggplot(Nd.long, aes(x = km, y = SNPsPerKb)) +
    geom_point(alpha = 0.25, shape = 1) +
    geom_smooth(method = "lm", linetype = 2, color = "black") +
    labs(x = "Geographic distance (km)", y = "Hamming distance (SNPs per Kb)", title = "b") +
    scale_y_continuous(breaks = c(1,4,7,10)) +
    annotate(
        geom = "text", 
        label = expr(paste("Mantel r = ", !!Nd.geo.mantel.stat, ", ", italic("P"), " = ", !!Nd.geo.mantel.sig, 1)),
        x = 1600,
        y = 0.75
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.06),
        axis.title.y = element_blank()
    )
p2

pdf("figures/pop_gen/IBD/IBD.pdf", width = 10, height = 3.5)
grid.arrange(p1,p2,ncol = 2)
dev.off()

png("figures/pop_gen/IBD/IBD.png", width = 10, height = 3.5, units = "in", res = 300)
grid.arrange(p1,p2,ncol = 2)
dev.off()

#duration
Nf.dur.mantel.stat = round(Nf.dur.mantel$statistic,2)
Nf.dur.mantel.sig = round(Nf.dur.mantel$signif,3)

p3 = ggplot(Nf.long, aes(x = durDif, y = SNPsPerKb)) +
    geom_point(alpha = 0.08, shape = 1) +
    geom_smooth(method = "lm", linetype = 1, color = "white", linewidth = 2) +
    geom_smooth(method = "lm", linetype = 1, color = "black") +
    labs(x = "Difference in infestation duration (years)", y = "Hamming distance (SNPs per Kb)", title = "a") +
    annotate(
        geom = "text", 
        label = expr(paste("Mantel r = ", !!Nf.dur.mantel.stat, ", ", italic("P"), " = ", !!Nf.dur.mantel.sig)),
        x = 61,
        y = 2.5
    ) +
    scale_y_continuous(breaks = c(2,3,4,5)) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.10)
    )
p3

Nd.dur.mantel.stat = round(Nd.dur.mantel$statistic,2)
Nd.dur.mantel.sig = round(Nd.dur.mantel$signif,2)


p4 = ggplot(Nd.long, aes(x = durDif, y = SNPsPerKb)) +
    geom_point(alpha = 0.25, shape = 1) +
    geom_smooth(method = "lm", linetype = 2, color = "black") +
    labs(x = "Difference in infestation duration (years)", y = "Hamming distance (SNPs per Kb)", title = "b") +
    scale_y_continuous(breaks = c(1,4,7,10)) +
    annotate(
        geom = "text", 
        label = expr(paste("Mantel r = ", !!Nd.dur.mantel.stat, ", ", italic("P"), " = ", !!Nd.dur.mantel.sig)),
        x = 69,
        y = 0.75
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.06),
        axis.title.y = element_blank()
    )
p4

pdf("figures/pop_gen/IBD/IB_durationInfection.pdf", width = 10, height = 3.5)
grid.arrange(p3,p4,ncol = 2)
dev.off()

png("figures/pop_gen/IBD/IB_durationInfection.png", width = 10, height = 3.5, units = "in", res = 300)
grid.arrange(p3,p4,ncol = 2)
dev.off()


####################
# four panel fig

p1 = ggplot(Nf.long, aes(x = km, y = SNPsPerKb)) +
    geom_point(alpha = 0.16, shape = 1) +
    geom_smooth(method = "lm", linetype = 1, color = "white", linewidth = 2) +
    geom_smooth(method = "lm", linetype = 1, color = "black") +
    scale_y_continuous(breaks = c(3,4,5)) +
    labs(x = "Geographic distance (km)", y = "Hamming distance (SNPs per Kb)", title = "a") +
    annotate(
        geom = "text", 
        label = expr(paste("Mantel r = ", !!Nf.geo.mantel.stat, ", ", italic("P"), " = ", !!Nf.geo.mantel.sig)),
        x = 1475,
        y = 2.5
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.10, vjust = -1),
        axis.title.x = element_blank()
    )
p1

p2 = ggplot(Nd.long, aes(x = km, y = SNPsPerKb)) +
    geom_point(alpha = 0.25, shape = 1) +
    geom_smooth(method = "lm", linetype = 2, color = "black") +
    labs(x = "Geographic distance (km)", y = "Hamming distance (SNPs per Kb)", title = "c") +
    scale_y_continuous(breaks = c(1,4,7,10)) +
    annotate(
        geom = "text", 
        label = expr(paste("Mantel r = ", !!Nd.geo.mantel.stat, ", ", italic("P"), " = ", !!Nd.geo.mantel.sig,1)),
        x = 1600,
        y = 0.75
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.115, vjust = -1)
    )


p3 = ggplot(Nf.long, aes(x = durDif, y = SNPsPerKb)) +
    geom_point(alpha = 0.16, shape = 1) +
    geom_smooth(method = "lm", linetype = 1, color = "white", linewidth = 2) +
    geom_smooth(method = "lm", linetype = 1, color = "black") +
    labs(x = "Difference in infestation duration (years)", y = "", title = "b") +
    annotate(
        geom = "text", 
        label = expr(paste("Mantel r = ", !!Nf.dur.mantel.stat, ", ", italic("P"), " = ", !!Nf.dur.mantel.sig)),
        x = 61,
        y = 2.5
    ) +
    scale_y_continuous(breaks = c(2,3,4,5)) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.04, vjust = -1),
        #axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank()
    )
p3


p4 = ggplot(Nd.long, aes(x = durDif, y = SNPsPerKb)) +
    geom_point(alpha = 0.25, shape = 1) +
    geom_smooth(method = "lm", linetype = 2, color = "black") +
    labs(x = "Difference in infestation duration (years)", y = "", title = "d") +
    scale_y_continuous(breaks = c(1,4,7,10)) +
    annotate(
        geom = "text", 
        label = expr(paste("Mantel r = ", !!Nd.dur.mantel.stat, ", ", italic("P"), " = ", !!Nd.dur.mantel.sig)),
        x = 69,
        y = 0.75
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.04, vjust = -1),
        #axis.title.y = element_blank(),
        axis.text.y = element_blank()
    )
p4


pdf("figures/pop_gen/IBD/IBD_durationInfection.four_panel.pdf", width = 10, height = 7)
grid.arrange(p1,p3,p2,p4,ncol = 2)
dev.off()
