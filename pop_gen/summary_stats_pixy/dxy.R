library(dplyr)
library(tidyr)
library(geosphere)
library(vegan)


sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
state_n = sample_metadata.Nf %>% group_by(state) %>% summarise(n = n())
state_n %>% print(n = Inf)


dxy = read.table("data/Nf/pixy/windowed_10kb/pixy_dxy.txt", header = T)
dxy_means = dxy %>%
    group_by(pop1, pop2) %>%
    summarize(dxy_mean = sum(count_diffs, na.rm = T)/sum(count_comparisons, na.rm = T))

dxy_means
#########
#########
#filter based on n
n_min = 4
low_n = state_n %>% filter(n < n_min) %>% pull(state)
dxy_means %>%
    filter(!pop1 %in% low_n & !pop2 %in% low_n) -> dxy_means.min_n

dxy_means.min_n

#############
dxy_mat = as.dist(xtabs(dxy_mean ~ pop2 + pop1, data.frame(dxy_means.min_n)))
head(dxy_mat)
dxy_mat
# the xtabs form does not work because of the ordering.
# we end up with zeros where there shouldn't be

dfr <- reshape(data.frame(dxy_means.min_n), direction="wide", idvar="pop1", timevar="pop2")
head(dfr)
dxy_means %>% filter(pop2 == "CT" & pop1 == "ME.N")
class(dfr)

dfr.ordered = dfr[order(rowSums(is.na(dfr)), decreasing = T),order(colSums(is.na(dfr)), decreasing = F)]
rownames(dfr.ordered) = dfr.ordered$pop1
#pop1 is the first col
dfr.ordered[,"pop1"] = NULL
#the colnames have extra text
colnames(dfr.ordered) = sub("dxy_mean.", "", colnames(dfr.ordered))
class(dfr.ordered)
dfr.ordered %>% lower.tri()
dfr.ordered %>% upper.tri()
dxy.mat = matrix(nrow = nrow(dfr.ordered)+1, ncol = ncol(dfr.ordered)+1)
dxy.mat[lower.tri(dxy.mat)] = dfr.ordered[!is.na(dfr.ordered)]

#full site order (1-16)
# need to add the last index of rownames to colnames
# (or the first of colnames to rownames) bc dist
sites_ordered = c(colnames(dfr.ordered), tail(rownames(dfr.ordered), 1))

colnames(dxy.mat) = sites_ordered
rownames(dxy.mat) = sites_ordered

dxy.dist = as.dist(dxy.mat)

###################
#setting up geo distance
site_dat = sample_metadata.Nf %>% 
    select(state, lat, lon) %>%
    filter(!state %in% low_n) %>% 
    distinct
row.names(site_dat) = site_dat$state
site_dat.ordered = site_dat[sites_ordered,]

Dgeo = distm(x = site_dat.ordered[,c("lon", "lat")], fun = distVincentyEllipsoid)
rownames(Dgeo) = site_dat.ordered$state
colnames(Dgeo) = site_dat.ordered$state
Dgeo = Dgeo/1000
as.dist(Dgeo)
#mantel
mantel(dxy.dist, as.dist(Dgeo))
#Mantel statistic r: 0.5948 
#Significance: 0.002 

#with 4 samples per site
#Mantel statistic r: 0.5872 
#Significance: 0.002 

#INCLUDING ALL SITES (EVEN WITH 1 REP)
#Mantel statistic r: 0.3724 
#      Significance: 0.007 


