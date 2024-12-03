library(dplyr)
library(vcfR)
#library(adegenet)
library(snpR)
# might need to import dartR, if the gl2snpR conversion fails try loading it


#metadata
sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")

#########
#prev used table of sample,state.name,lat,lon
ind.metrics = sample_metadata.Nf %>% select(Sequence_label, state, lat, lon)
colnames(ind.metrics) = c("sampID", "pop", "lat", "lon")
sub("\\.", "_", ind.metrics$pop) -> ind.metrics$pop #snpR does not accept "." ... beyond annoying

row.names(ind.metrics) = ind.metrics$sampID

#filtered VCF
# bialleles only
# we are now also excluding singletons bc this was causing problems with the conversion to snpR
vcf <- read.vcfR("data/Nf/final_tables/rm_dups/FINAL_snp.biallele.vcf.gz", verbose = FALSE)

gt = extract.gt(vcf, element='GT', as.numeric=TRUE)
gt.pos_list = strsplit(row.names(gt), "_")
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

gt = data.frame(gt)
nrow(gt)
nrow(gt.pos_df)
gt[is.na(gt)] = "NN"
head(gt)


dat = import.snpR.data(gt, snp.meta = gt.pos_df, sample.meta = ind.metrics)
dat

rm(vcf)
gc()


#bi-allelic snpRdata with 424811 SNPs and 115 samples.
#Calculated statistics can be accessed via get.snpR.stats()
#
# note that this is about half of the SNPs we get when including mac==1
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

#foo = get.snpR.stats(dat)
str(foo)
foo@weighted.means

tajD = foo@weighted.means %>% 
    filter(facet == "pop" & !is.na(global_D)) %>%
    select(subfacet, global_ws.theta, global_ts.theta, global_D, global_num_seg)
colnames(tajD)[1] = "state"
tajD$state = sub("_", ".", tajD$state)
write.csv(tajD, "data/Nf/pixy/tajD.snpR.csv", row.names = F, quote = F)


