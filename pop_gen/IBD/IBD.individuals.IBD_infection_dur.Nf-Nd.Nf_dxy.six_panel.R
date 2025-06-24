library(dplyr)
library(ggplot2)
library(gridExtra)
library(vegan)
library(geosphere)
library(rlang)
library(gtable)
source("library/ggplot_theme.txt")

Dgen.Nf = read.table("data/Nf/final_tables/rm_dups/FINAL_invariant.IBD_analyses.dist", header = F) 
dist.ID.Nf = read.table("data/Nf/final_tables/rm_dups/FINAL_invariant.IBD_analyses.dist.id", header = F)
rownames(Dgen.Nf) = dist.ID.Nf[,1]
colnames(Dgen.Nf) = dist.ID.Nf[,2]
sum(is.na(Dgen.Nf))
Dgen.Nf = Dgen.Nf %>% as.dist()
Dgen.Nf
sum(is.na(Dgen.Nf))

Dgen.Nd = read.table("data/Nd/final_tables/rm_dups/FINAL_invariant.IBD_analyses.dist", header = F) 
dist.ID.Nd = read.table("data/Nd/final_tables/rm_dups/FINAL_invariant.IBD_analyses.dist.id", header = F)
rownames(Dgen.Nd) = dist.ID.Nd[,1]
colnames(Dgen.Nd) = dist.ID.Nd[,2]
Dgen.Nd = Dgen.Nd %>% as.dist()
head(Dgen.Nd)
sum(is.na(Dgen.Nd))

#convert to nt difs per kb
Dgen.Nf = Dgen.Nf / (41040857 / 1000)
sum(is.na(Dgen.Nf))
Dgen.Nd = Dgen.Nd / (38535154 / 1000)
sum(is.na(Dgen.Nd))

#get order of samples
dist.order.Nf = as.matrix(Dgen.Nf ) %>% rownames
dist.order.Nd = as.matrix(Dgen.Nd ) %>% rownames 

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

#dur.Nf$lat[dur.Nf$collection_period == "early"] = NA
#dur.Nf$lon[dur.Nf$collection_period == "early"] = NA
#dur.Nf$dur_inf[dur.Nf$collection_period == "early"] = NA

#dur.Nd$lat[dur.Nd$collection_period == "early"] = NA
#dur.Nd$lon[dur.Nd$collection_period == "early"] = NA
#dur.Nd$dur_inf[dur.Nd$collection_period == "early"] = NA

#calculate geographic distance and age dif
Nf.Dgeo <- distm(x = dur.Nf[,c("lon", "lat")], fun = distVincentyEllipsoid)
Nd.Dgeo <- distm(x = dur.Nd[,c("lon", "lat")], fun = distVincentyEllipsoid)
Nf.Ddur <- dist(x = dur.Nf$dur_inf) %>% as.matrix
Nd.Ddur <- dist(x = dur.Nd$dur_inf) %>% as.matrix
#set names
rownames(Nf.Ddur) = dur.Nf$samp
colnames(Nf.Ddur) = dur.Nf$samp
rownames(Nf.Dgeo) = dur.Nf$samp
colnames(Nf.Dgeo) = dur.Nf$samp

rownames(Nd.Ddur) = dur.Nd$samp
colnames(Nd.Ddur) = dur.Nd$samp
rownames(Nd.Dgeo) = dur.Nd$samp
colnames(Nd.Dgeo) = dur.Nd$samp

#converting geo to dist class now before setting within site comps to NA
Nf.Dgeo = as.dist(Nf.Dgeo/1000)
Nd.Dgeo = as.dist(Nd.Dgeo/1000)

#set zero dists to NA (within site comps) AND selfs (diag)
sum(Nf.Dgeo == 0, na.rm = T)
#341
#no change with early
Nf.Dgeo[Nf.Dgeo == 0] = NA
sum(Nf.Dgeo == 0, na.rm = T)

sum(Nd.Dgeo == 0, na.rm = T)
#20
# 23 with early
Nd.Dgeo[Nd.Dgeo == 0] = NA

#make dur dist class
Nf.Ddur = as.dist(Nf.Ddur)
Nd.Ddur = as.dist(Nd.Ddur)

sum(is.na(Dgen.Nf))
sum(is.na(Nf.Ddur))
sum(is.na(Nf.Dgeo))
sum(is.na(Dgen.Nd))
sum(is.na(Nd.Ddur))
sum(is.na(Nd.Dgeo))

# match NAs in gen dists and durs
Dgen.Nf[is.na(Nf.Dgeo)] = NA
Nf.Ddur[is.na(Nf.Dgeo)] = NA
Dgen.Nd[is.na(Nd.Dgeo)] = NA
Nd.Ddur[is.na(Nd.Dgeo)] = NA

sum(is.na(Dgen.Nf))
sum(is.na(Nf.Ddur))
sum(is.na(Nf.Dgeo))
sum(is.na(Dgen.Nd))
sum(is.na(Nd.Ddur))
sum(is.na(Nd.Dgeo))

####################
#Mantel geo dist
mantel(Dgen.Nf, Nf.Dgeo, na.rm = T)
#Mantel statistic r: 0.2467
#Significance: 0.001
#with early
#Mantel statistic r: 0.2538 
#      Significance: 0.001 

mantel(Dgen.Nd, Nd.Dgeo, na.rm = T)
#Mantel statistic r: -0.09242
#Significance: 0.788
#with early
#Mantel statistic r: -0.1205 
#      Significance: 0.904

mantel(Dgen.Nf, log(Nf.Dgeo), na.rm = T)
#Mantel statistic r: 0.2304
#Significance: 0.001
mantel(Dgen.Nd, log(Nd.Dgeo), na.rm = T)
#Mantel statistic r: -0.1145
#Significance: 0.887

####################
#Mantel duration dist
mantel(Dgen.Nf, Nf.Ddur, na.rm = T)
#Mantel statistic r: 0.2144 
#Significance: 0.001 
##with early
#Mantel statistic r: 0.2045 
#      Significance: 0.001 
mantel(Dgen.Nd, Nd.Ddur, na.rm = T)
#Mantel statistic r: 0.01163 
#Significance: 0.454 
##with early
#Mantel statistic r: 0.02946 
#      Significance: 0.337 

#########################
#long format for plotting
#set all diags and upper.tris to NA before long format

#geo
Nf.Dgeo.mat = Nf.Dgeo %>% as.matrix
diag(Nf.Dgeo.mat) = NA
Nf.Dgeo.mat[upper.tri(Nf.Dgeo.mat)] = NA
Nf.Dgeo.long = reshape2::melt(Nf.Dgeo.mat)
Nf.Dgeo.long %>% filter(Var1 == "NG2" & Var2 == "NG1") %>% nrow

Nd.Dgeo.mat = Nd.Dgeo %>% as.matrix
diag(Nd.Dgeo.mat) = NA
Nd.Dgeo.mat[upper.tri(Nd.Dgeo.mat)] = NA
Nd.Dgeo.long = reshape2::melt(Nd.Dgeo.mat)

#gen
Nf.Dgen.mat = Dgen.Nf %>% as.matrix
diag(Nf.Dgen.mat) = NA
Nf.Dgen.mat[upper.tri(Nf.Dgen.mat)] = NA
Nf.Dgen.long = reshape2::melt(Nf.Dgen.mat)
#checks
sum(is.na(Dgen.Nf))
sum(is.na(Nf.Dgen.mat))
sum(is.na(Nf.Dgen.long$value))

Nd.Dgen.mat = Dgen.Nd %>% as.matrix
diag(Nd.Dgen.mat) = NA
Nd.Dgen.mat[upper.tri(Nd.Dgen.mat)] = NA
Nd.Dgen.long = reshape2::melt(Nd.Dgen.mat)

#duration
Nf.Ddur.mat = Nf.Ddur %>% as.matrix
diag(Nf.Ddur.mat) = NA
Nf.Ddur.mat[upper.tri(Nf.Ddur.mat)] = NA
Nf.Ddur.long = reshape2::melt(Nf.Ddur.mat)
#checks
Nf.Ddur.long %>% filter(Var1 == "NG121" & Var2 == "NG4") 
Nf.Ddur.long %>% filter(Var1 == "NG4" & Var2 == "NG121")
Nf.Dgeo.long %>% filter(Var1 == "NG121" & Var2 == "NG4") 
Nf.Dgeo.long %>% filter(Var1 == "NG4" & Var2 == "NG121")
Nf.Dgen.long %>% filter(Var1 == "NG121" & Var2 == "NG4")
Nf.Dgen.long %>% filter(Var1 == "NG4" & Var2 == "NG121")

Nd.Ddur.mat = Nd.Ddur %>% as.matrix
diag(Nd.Ddur.mat) = NA
Nd.Ddur.mat[upper.tri(Nd.Ddur.mat)] = NA
Nd.Ddur.long = reshape2::melt(Nd.Ddur.mat %>% as.matrix)

#check all equal NAs
sum(is.na(Nf.Dgen.long$value)) == sum(is.na(Nf.Ddur.long$value)) 
sum(is.na(Nf.Dgen.long$value)) == sum(is.na(Nf.Dgeo.long$value))

sum(is.na(Nd.Dgen.long$value)) == sum(is.na(Nd.Ddur.long$value)) 
sum(is.na(Nd.Dgen.long$value)) == sum(is.na(Nd.Dgeo.long$value))

#name
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
nrow(Nd.long)
nrow(Nd.Dgeo.long)
nrow(Nd.Dgen.long)

#rm NAs
Nf.long = Nf.long[!is.na(Nf.long$km),]
Nd.long = Nd.long[!is.na(Nd.long$km),]

###################################
###################################
###################################
# DXY versus geo dist 
#
dxy.dist = readRDS("data/Nf/pixy/windowed_10kb/dxy_dist.rds")
dxy.mat = as.matrix(dxy.dist)
dxy.mat
nrow(dxy.mat)
ncol(dxy.mat)

state_n = sample_metadata.Nf %>% group_by(state) %>% summarise(n = n())
state_n %>% print(n = Inf)
n_min = 4
low_n = state_n %>% filter(n < n_min) %>% pull(state)

###################
#setting up geo distance
site_dat = sample_metadata.Nf %>% 
    select(state, lat, lon) %>%
    filter(!state %in% low_n) %>% 
    distinct
nrow(site_dat) == nrow(dxy.mat)
row.names(site_dat) = site_dat$state
site_dat.ordered = site_dat[rownames(dxy.mat),]

Dgeo.NfPop = distm(x = site_dat.ordered[,c("lon", "lat")], fun = distVincentyEllipsoid)
rownames(Dgeo.NfPop) = site_dat.ordered$state
colnames(Dgeo.NfPop) = site_dat.ordered$state
Dgeo.NfPop = Dgeo.NfPop/1000
as.dist(Dgeo.NfPop)
#mantel
Nf.geoPop.mantel = mantel(dxy.dist, as.dist(Dgeo.NfPop))
Nf.geoPop.mantel
#Mantel statistic r: 0.5872 
#Significance: 0.002 
mantel(dxy.dist, log(as.dist(Dgeo.NfPop)))
#Mantel statistic r: 0.5574 
#Significance: 0.001 

###################
#setting up dur distance
#
site_dur = sample_metadata.Nf %>% 
    select(state, duration_infection) %>%
    filter(!state %in% low_n) %>% 
    distinct
nrow(site_dur) == nrow(dxy.mat)
row.names(site_dur) = site_dur$state
site_dur.ordered = site_dur[rownames(dxy.mat),]

Nf.durPop <- dist(x = site_dur.ordered$duration_infection) %>% as.matrix
rownames(Nf.durPop) = site_dur.ordered$state
colnames(Nf.durPop) = site_dur.ordered$state

Nf.durPop.mantel = mantel(dxy.dist, as.dist(Nf.durPop))
Nf.durPop.mantel
#Mantel statistic r: 0.5729 
#Significance: 0.001 

#long formats
diag(dxy.mat) = NA
dxy.mat[upper.tri(dxy.mat)] = NA
Nf.distPop.long = reshape2::melt(dxy.mat)
sum(is.na(Nf.distPop.long$value))
sum(!is.na(Nf.distPop.long$value))
Nf.distPop.long = Nf.distPop.long[!is.na(Nf.distPop.long$value),]
colnames(Nf.distPop.long) = c("pop1", "pop2", "dxy")

diag(Dgeo.NfPop) = NA
Dgeo.NfPop[upper.tri(Dgeo.NfPop)] = NA
Nf.geoPop.long = reshape2::melt(Dgeo.NfPop)
sum(is.na(Nf.geoPop.long$value))
sum(!is.na(Nf.geoPop.long$value))
Nf.geoPop.long = Nf.geoPop.long[!is.na(Nf.geoPop.long$value),]
colnames(Nf.geoPop.long) = c("pop1", "pop2", "km")

diag(Nf.durPop) = NA
Nf.durPop[upper.tri(Nf.durPop)] = NA
Nf.durPop.long = reshape2::melt(Nf.durPop %>% as.matrix)
sum(is.na(Nf.durPop.long$value))
sum(!is.na(Nf.durPop.long$value))
Nf.durPop.long = Nf.durPop.long[!is.na(Nf.durPop.long$value),]
colnames(Nf.durPop.long) = c("pop1", "pop2", "dur_dist")

Nf_pop.long = left_join(Nf.distPop.long, Nf.geoPop.long, by = c("pop1", "pop2")) %>%
    left_join(., Nf.durPop.long, by = c("pop1", "pop2"))

Nf_pop.geo.mantel.stat = round(Nf.geoPop.mantel$statistic,2)
Nf_pop.geo.mantel.sig = round(Nf.geoPop.mantel$signif,3)
Nf_pop.dur.mantel.stat = round(Nf.durPop.mantel$statistic,2)
Nf_pop.dur.mantel.sig = round(Nf.durPop.mantel$signif,3)

######################
#plot
#
#stats from indv comps
Nf.geo.mantel = mantel(Dgen.Nf, Nf.Dgeo, na.rm = T)
Nd.geo.mantel = mantel(Dgen.Nd, Nd.Dgeo, na.rm = T)
Nf.dur.mantel = mantel(Dgen.Nf, Nf.Ddur, na.rm = T)
Nd.dur.mantel = mantel(Dgen.Nd, Nd.Ddur, na.rm = T)

Nf.geo.mantel.stat = round(Nf.geo.mantel$statistic,2)
Nf.geo.mantel.sig = round(Nf.geo.mantel$signif,3)
Nd.geo.mantel.stat = round(Nd.geo.mantel$statistic,2)
Nd.geo.mantel.sig = round(Nd.geo.mantel$signif,2)

#duration
Nf.dur.mantel.stat = round(Nf.dur.mantel$statistic,2)
Nf.dur.mantel.sig = round(Nf.dur.mantel$signif,3)
Nd.dur.mantel.stat = round(Nd.dur.mantel$statistic,2)
Nd.dur.mantel.sig = round(Nd.dur.mantel$signif,2)


####################
# PLOTTING
####################
####################
p1 = ggplot(Nf.long, aes(x = km, y = SNPsPerKb)) +
    geom_point(alpha = 0.16, shape = 1) +
    geom_smooth(method = "lm", linetype = 1, color = "white", linewidth = 2) +
    geom_smooth(method = "lm", linetype = 1, color = "black") +
    scale_y_continuous(breaks = c(3,4,5)) +
    labs(
        x = "Geographic distance (km)", 
        y = expression(paste("Hamming distance (SNPs Kb"^-1,")")), 
        title = "a"
    ) +
    annotate(
        geom = "text", 
        label = expr(paste("Mantel r = ", !!Nf.geo.mantel.stat, ", ", italic("P"), " = ", !!Nf.geo.mantel.sig)),
        x = max(Nf.long$km),#1475,
        y = min(Nf.long$SNPsPerKb),#2.5,
        hjust = 1
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.10, vjust = -1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 1.5, b = 0, l = 0))
    )
p1

p2 = ggplot(Nd.long, aes(x = km, y = SNPsPerKb)) +
    geom_point(alpha = 0.35, shape = 1) +
    geom_smooth(method = "lm", linetype = 2, color = "black") +
    labs(
        x = "Geographic distance (km)", 
        y = expression(paste("Hamming distance (SNPs Kb"^-1,")")), 
        title = "c"
    ) +
    scale_y_continuous(breaks = c(1,4,7,10)) +
    annotate(
        geom = "text", 
        label = expr(paste("Mantel r = ", !!Nd.geo.mantel.stat, ", ", italic("P"), " = ", !!Nd.geo.mantel.sig)),
        x = max(Nd.long$km),#1600,
        y = min(Nd.long$SNPsPerKb),#0.75,
        hjust = 1
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.10, vjust = -1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 1.5, b = 0, l = 0))
    )
p2

p3 = ggplot(Nf.long, aes(x = durDif, y = SNPsPerKb)) +
    geom_point(alpha = 0.16, shape = 1) +
    geom_smooth(method = "lm", linetype = 1, color = "white", linewidth = 2) +
    geom_smooth(method = "lm", linetype = 1, color = "black") +
    labs(x = "Difference in infestation duration (years)", y = "", title = "b") +
    annotate(
        geom = "text", 
        label = expr(paste("Mantel r = ", !!Nf.dur.mantel.stat, "0, ", italic("P"), " = ", !!Nf.dur.mantel.sig)),
        x = max(Nf.long$durDif),#61,
        y = min(Nf.long$SNPsPerKb),#2.5,
        hjust = 1
    ) +
    scale_y_continuous(breaks = c(2,3,4,5)) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.03, vjust = -1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(l = 12, r = 5.5, b = 5.5, t = 5.5)
    )
p3


p4 = ggplot(Nd.long, aes(x = durDif, y = SNPsPerKb)) +
    geom_point(alpha = 0.35, shape = 1) +
    geom_smooth(method = "lm", linetype = 2, color = "black") +
    labs(x = "Difference in infestation duration (years)", y = "", title = "d") +
    scale_y_continuous(breaks = c(1,4,7,10)) +
    annotate(
        geom = "text", 
        label = expr(paste("Mantel r = ", !!Nd.dur.mantel.stat, ", ", italic("P"), " = ", !!Nd.dur.mantel.sig)),
        x = max(Nd.long$durDif),#69,
        y = min(Nd.long$SNPsPerKb),#0.75,
        hjust = 1
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.03, vjust = -1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(l = 12, r = 5.5, b = 5.5, t = 5.5)
    )
p4

p5 = ggplot(Nf_pop.long, aes(x = km, y = dxy*10^3)) +
    geom_smooth(method = "lm", linetype = 1, color = "black", se = T) +
    geom_point(shape = 1) +
    labs(x = "Geographic distance (km)", 
         y = expression(paste(d[xy]%*%10^3, " (SNPs Kb"^-1,")")), 
         title = "e"
    ) +
    annotate(
        geom = "text", 
        label = expr(paste("Mantel r = ", !!Nf_pop.geo.mantel.stat, ", ", italic("P"), " = ", !!Nf_pop.geo.mantel.sig)),
        x = max(Nf_pop.long$km),#1410,
        y = min(Nf_pop.long$dxy*10^3),
        hjust = 1
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.10, vjust = -1),
        axis.title.y = element_text(margin = margin(t = 0, r = 1.5, b = 0, l = 0))
    )
p5    

p6 = ggplot(Nf_pop.long, aes(x = dur_dist, y = dxy*10^3)) +
    geom_smooth(method = "lm", linetype = 1, color = "black", se = T) +
    geom_point(shape = 1) +
    labs(
        x = "Difference in infestation duration (years)", 
         y = "", 
         title = "f"
    ) +
    annotate(
        geom = "text", 
        label = expr(paste("Mantel r = ", !!Nf_pop.dur.mantel.stat, ", ", italic("P"), " = ", !!Nf_pop.dur.mantel.sig)),
        x = max(Nf_pop.long$dur_dist),#50
        y = min(Nf_pop.long$dxy*10^3),
        hjust = 1
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.03, vjust = -1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(l = 12, r = 5.5, b = 5.5, t = 5.5)
    )
p6    

gp1 = ggplotGrob(p1)
gp2 = ggplotGrob(p2)
gp3 = ggplotGrob(p3)
gp4 = ggplotGrob(p4)
gp5 = ggplotGrob(p5)
gp6 = ggplotGrob(p6)

gcol1 = rbind(gp1, gp2, gp5)
gcol2 = rbind(gp3, gp4, gp6)

pdf("figures/pop_gen/IBD/IBD_durationInfection.six_panel.include_early.pdf", width = 10, height = 10)
grid.arrange(gcol1, gcol2, ncol = 2, widths = c(0.525, 0.475))
dev.off()

png("figures/pop_gen/IBD/IBD_durationInfection.six_panel.include_early.png", width = 10, height = 10, units = "in", res = 300)
grid.arrange(gcol1, gcol2, ncol = 2, widths = c(0.525, 0.475))
dev.off()
