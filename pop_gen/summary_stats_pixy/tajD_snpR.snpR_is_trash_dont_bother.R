library(dplyr)
library(vcfR)
#library(adegenet)
library(snpR)
# might need to import dartR, if the gl2snpR conversion fails try loading it


#metadata
sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
state_n = sample_metadata.Nf %>% group_by(state) %>% summarise(n = n())
state_n %>% print(n = Inf)
#filter based on n
n_min = 4
low_n = state_n %>% filter(n < n_min) %>% pull(state)

#########
#prev used table of sample,state.name,lat,lon
ind.metrics = sample_metadata.Nf %>% 
    select(Sequence_label, state, lat, lon) %>% 
    filter(!state %in% low_n)
colnames(ind.metrics) = c("sampID", "pop", "lat", "lon")
sub("\\.", "_", ind.metrics$pop) -> ind.metrics$pop #snpR does not accept "." ... beyond annoying

row.names(ind.metrics) = ind.metrics$sampID

#filtered VCF
# bialleles only
# 
vcf <- read.vcfR("data/Nf/final_tables/rm_dups/FINAL_snp.biallele.vcf.gz", verbose = FALSE)

gt = extract.gt(vcf, element='GT', as.numeric=TRUE)
gt.low_n = gt[,colnames(gt) %in% ind.metrics$sampID]
ncol(gt.low_n)
gt.pos_list = strsplit(row.names(gt.low_n), "_")
gt.pos_df = data.frame(
    chr = paste(
        lapply(gt.pos_list, function(x) x[1]) %>% unlist,
        lapply(gt.pos_list, function(x) x[2]) %>% unlist,
        sep = "_"
    ),
    position = lapply(gt.pos_list, function(x) x[3]) %>% unlist,
    stringsAsFactors = F
) #note these column names must be chr and position for snpR
#because theres no argument to indicate chr and pos in the tajimaD func

gt = data.frame(gt.low_n)
nrow(gt)
nrow(gt.pos_df)
gt[is.na(gt)] = "NN"
head(gt)


dat = import.snpR.data(gt, snp.meta = gt.pos_df, sample.meta = ind.metrics)
dat

rm(vcf)
gc()


#bi-allelic snpRdata with 969175 SNPs and 115 samples.
#Calculated statistics can be accessed via get.snpR.stats()
#
# 
# 
foo = calc_tajimas_d(
    dat,
    facets = c("pop", "chr"),
    #sigma = 100, #sliding window size in kb
    #step = 200, #defulat = sigma*2 (non-overlapping windows)
    par = 5, #number cores
    #triple_sigma = F, #this is for smoothing, calcualtes average in window of 3xsigma 
    global = T #use this to just calculate global
)

foo_seg = calc_seg_sites(
    dat, 
    facets = c("pop", "chr"),
    rarefaction = T
)
        
                         
#foo = get.snpR.stats(dat)
str(foo)
foo@weighted.means

tajD = foo@weighted.means %>% 
    filter(facet == "pop" & !is.na(global_D)) %>%
    select(subfacet, global_ws.theta, global_ts.theta, global_D, global_num_seg)
colnames(tajD)[1] = "state"
tajD$state = sub("_", ".", tajD$state)
write.csv(tajD, "data/Nf/pixy/tajD.snpR.csv", row.names = F, quote = F)


