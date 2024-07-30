library(pegas)
library(dplyr)
library(vcfR)
library(adegenet)
library(poppr)

dist.raw.Nf = readRDS("data/Nf/IBD/hamming_dist.rds")
dist.raw.Nd = readRDS("data/Nd/IBD/hamming_dist.rds")

dist.order.Nf = as.matrix(dist.raw.Nf ) %>% rownames 
dist.order.Nd = as.matrix(dist.raw.Nd ) %>% rownames 

#metadata
sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
sample_metadata.Nd = read.csv("data/sample_metadata/Nd_filtered.lat_lon_dur_inf.csv")
rownames(sample_metadata.Nf) = sample_metadata.Nf$Sequence_label
rownames(sample_metadata.Nd) = sample_metadata.Nd$Sequence_label

hier.Nf = data.frame(
    samp = sample_metadata.Nf[dist.order.Nf, "Sequence_label"],
    site = sample_metadata.Nf[dist.order.Nf, "state"],
    stringsAsFactors = T
)
rownames(hier.Nf) = hier.Nf$samp

hier.Nd = data.frame(
    samp = sample_metadata.Nd[dist.order.Nd, "Sequence_label"],
    site = sample_metadata.Nd[dist.order.Nd, "state"],
    stringsAsFactors = T
)
rownames(hier.Nd) = hier.Nd$samp

set.seed(12345)

amova.pegas.Nf = pegas::amova(
    dist.raw.Nf ~ site, 
    is.squared = F, 
    data = hier.Nf, 
    nperm = 999 #set this to 999
)

amova.pegas.Nd = pegas::amova(
    dist.raw.Nd ~ site, 
    is.squared = F, 
    data = hier.Nd, 
    nperm = 999 #set this to 999
)


#################################
#################################
# To plot w/n vb/n dists; convert
# dist obj to a matrix, classify 
# each cell as within or b/n,
# then plot in a histogram