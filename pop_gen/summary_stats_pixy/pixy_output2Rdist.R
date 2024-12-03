library(dplyr)
library(tidyr)
library(geosphere)
library(vegan)


sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
state_n = sample_metadata.Nf %>% group_by(state) %>% summarise(n = n())
state_n %>% print(n = Inf)
#filter based on n
n_min = 3
low_n = state_n %>% filter(n < n_min) %>% pull(state)

########################################################
########################################################
#first process dxy
dxy = read.table("data/Nf/pixy/pixy_dxy.txt", header = T)
dxy_means = dxy %>%
    group_by(pop1, pop2) %>%
    summarize(dxy_mean = sum(count_diffs, na.rm = T)/sum(count_comparisons, na.rm = T))

dxy_means
#########
#########
dxy_means %>%
    filter(!pop1 %in% low_n & !pop2 %in% low_n) -> dxy_means.min_n

dxy_means.min_n

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
dxy.mat[lower.tri(dxy.mat)] %>% length
dfr.ordered[!is.na(dfr.ordered)] %>% length
dxy.mat[lower.tri(dxy.mat)] = dfr.ordered[!is.na(dfr.ordered)]

#full site order (1-16)
# need to add the last index of rownames to colnames
# (or the first of colnames to rownames) bc pixy outputs *unique* pairwise comps 
# i.e., a lower/upper triangle comps matrix
sites_ordered = c(colnames(dfr.ordered), tail(rownames(dfr.ordered), 1))

colnames(dxy.mat) = sites_ordered
rownames(dxy.mat) = sites_ordered

dxy.dist = as.dist(dxy.mat)
saveRDS(object = dxy.dist, file = "data/Nf/pixy/dxy_dist.rds")

########################################################
########################################################
#process fst

fst = read.table("data/Nf/pixy/pixy_fst.txt", header = T)
fst_means = fst %>%
    group_by(pop1, pop2) %>%
    summarize(fst_mean = sum(count_diffs, na.rm = T)/sum(count_comparisons, na.rm = T))

fst_means
#########
#########
fst_means %>%
    filter(!pop1 %in% low_n & !pop2 %in% low_n) -> fst_means.min_n

fst_means.min_n

dfr <- reshape(data.frame(fst_means.min_n), direction="wide", idvar="pop1", timevar="pop2")
head(dfr)
fst_means %>% filter(pop2 == "CT" & pop1 == "ME.N")
class(dfr)

dfr.ordered = dfr[order(rowSums(is.na(dfr)), decreasing = T),order(colSums(is.na(dfr)), decreasing = F)]
rownames(dfr.ordered) = dfr.ordered$pop1
#pop1 is the first col
dfr.ordered[,"pop1"] = NULL
#the colnames have extra text
colnames(dfr.ordered) = sub("fst_mean.", "", colnames(dfr.ordered))
class(dfr.ordered)
dfr.ordered %>% lower.tri()
dfr.ordered %>% upper.tri()
fst.mat = matrix(nrow = nrow(dfr.ordered)+1, ncol = ncol(dfr.ordered)+1)
fst.mat[lower.tri(fst.mat)] %>% length
dfr.ordered[!is.na(dfr.ordered)] %>% length
fst.mat[lower.tri(fst.mat)] = dfr.ordered[!is.na(dfr.ordered)]

#full site order (1-16)
# need to add the last index of rownames to colnames
# (or the first of colnames to rownames) bc pixy outputs *unique* pairwise comps 
# i.e., a lower/upper triangle comps matrix
sites_ordered = c(colnames(dfr.ordered), tail(rownames(dfr.ordered), 1))

colnames(fst.mat) = sites_ordered
rownames(fst.mat) = sites_ordered

fst.dist = as.dist(fst.mat)
saveRDS(object = fst.dist, file = "data/Nf/pixy/fst_dist.rds")

