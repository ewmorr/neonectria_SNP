library(dplyr)
library(tidyr)
library(geosphere)
library(vegan)


########################################################
########################################################
#first process Nf
dxy = read.table("data/Nf/pixy/contig_all_dif_pop/pixy_dxy.txt", header = T)
dxy_means = dxy %>%
    group_by(pop1, pop2) %>%
    summarize(dxy_mean = sum(count_diffs, na.rm = T)/sum(count_comparisons, na.rm = T))

dxy_means
#########
#########
dfr <- reshape(data.frame(dxy_means), direction="wide", idvar="pop1", timevar="pop2")
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
saveRDS(object = dxy.dist, file = "data/Nf/pixy/contig_all_dif_pop/dxy_dist.rds")

########################################################
########################################################
#first process Nd
dxy = read.table("data/Nd/pixy/contig_all_dif_pop/pixy_dxy.txt", header = T)
dxy_means = dxy %>%
    group_by(pop1, pop2) %>%
    summarize(dxy_mean = sum(count_diffs, na.rm = T)/sum(count_comparisons, na.rm = T))

dxy_means
#########
#########
dfr <- reshape(data.frame(dxy_means), direction="wide", idvar="pop1", timevar="pop2")
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
saveRDS(object = dxy.dist, file = "data/Nd/pixy/contig_all_dif_pop/dxy_dist.rds")

########################################################
########################################################
#first process Nc
dxy = read.table("data/Nc/pixy/contig_all_dif_pop/pixy_dxy.txt", header = T)
dxy_means = dxy %>%
    group_by(pop1, pop2) %>%
    summarize(dxy_mean = sum(count_diffs, na.rm = T)/sum(count_comparisons, na.rm = T))

dxy_means
#########
#########
dfr <- reshape(data.frame(dxy_means), direction="wide", idvar="pop1", timevar="pop2")
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
saveRDS(object = dxy.dist, file = "data/Nc/pixy/contig_all_dif_pop/dxy_dist.rds")

