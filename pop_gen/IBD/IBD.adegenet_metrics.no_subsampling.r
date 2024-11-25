set.seed(12345)
library(vcfR)
library(dplyr)
library(adegenet)
library(geosphere)
library(vegan)

#premise conda env
#mamba create --name R-pop_gen conda-forge::r-tidyr conda-forge::r-dplyr bioconda::r-adegenet bioconda::r-vcfr
#conda activate R-pop_gen

#metadata
#sample_metadata.Nf = read.csv("~/Nf_pop_IBD_11182024/Nf_filtered.lat_lon_dur_inf.csv")
sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")

#########
#prev used table of sample,state.name,lat,lon
ind.metrics = sample_metadata.Nf %>% select(Sequence_label, state, lat, lon)

#filtered VCF
#vcf <- read.vcfR("~/Nf_pop_IBD_11182024/FINAL_snp.IBD_analyses.vcf.gz", verbose = FALSE)
vcf <- read.vcfR("data/Nf/final_tables/rm_dups/FINAL_snp.IBD_analyses.vcf.gz", verbose = FALSE)
gl = vcfR2genlight(vcf)
#In vcfR2genlight(vcf) : Found 47046 loci with more than two alleles.
#Objects of class genlight only support loci with two alleles.
#47046 loci will be omitted from the genlight object.
rm(vcf)
gc()


row.names(ind.metrics) = ind.metrics$Sequence_label

#Set metadata in genLight
gl@other$ind.metrics = ind.metrics[gl@ind.names,]
gl@other$latlong = ind.metrics[,3:4]
gl@other$latlong
gl@pop = as.factor(ind.metrics$state) #need to set pop for the ibd test to work
gl@ploidy = rep(as.integer(1), nInd(gl))

#SOME SITES HAVE SMALL SAMPLE SIZE
#Set min sample size to 3
low_n_states = ( (gl@other$ind.metrics %>% group_by(state) %>% summarize(n = n()) ) %>% filter(n < 3) )$state

keep.ind.list = data.frame(gl@ind.names, pop.gl = pop(gl)) %>% filter(!pop.gl %in% c(low_n_states))


#start.time <- Sys.time()
y = as.data.frame(gl)
gi = df2genind(y+1, ploidy=1)
gi@pop = gl@pop

gi.subset = gi[keep.ind.list$gl.ind.names]
to_remove = names(gi.subset@loc.n.all[gi.subset@loc.n.all == 1 ])
rm_indx = which(names(gi.subset@loc.n.all) %in% to_remove)
gi.subset.rm = gi.subset[loc=-rm_indx]
gp = genind2genpop(gi.subset.rm)

Dgen.2 <- dist.genpop(gp,method=2)

saveRDS(Dgen.2, "data/Nf/IBD/Nf.DSCE.no_subsample.rds")

#setting up distance
site_dat = ind.metrics %>% 
    filter(!state %in% low_n_states) %>% 
    select(!Sequence_label) %>%
    distinct
row.names(site_dat) = site_dat$state
site_dat.ordered = site_dat[row.names(gp@tab),]

Dgeo = distm(x = site_dat.ordered[,c("lon", "lat")], fun = distVincentyEllipsoid)
rownames(Dgeo) = site_dat.ordered$state
colnames(Dgeo) = site_dat.ordered$state
Dgeo = Dgeo/1000

#mantel
mantel(Dgen.2, Dgeo)

#mantel no VA
##############
#need to fully recalc without VA
# not just remove from the matrices
# 

keep.ind.list.no_VA = data.frame(gl@ind.names, pop.gl = pop(gl)) %>% filter(!pop.gl %in% c(low_n_states, "VA"))

gi.subset.no_VA = gi.subset[keep.ind.list.no_VA$gl.ind.names]
to_remove = names(gi.subset.no_VA@loc.n.all[gi.subset.no_VA@loc.n.all == 1 ])
rm_indx = which(names(gi.subset.no_VA@loc.n.all) %in% to_remove)
gi.subset.rm = gi.subset.no_VA[loc=-rm_indx]
gp = genind2genpop(gi.subset.rm)

Dgen.2 <- dist.genpop(gp,method=2)



site_dat.no_VA = ind.metrics %>% 
    filter(!state %in% c(low_n_states, "VA")) %>% 
    select(!Sequence_label) %>%
    distinct
row.names(site_dat.no_VA) = site_dat.no_VA$state
site_dat.ordered = site_dat[row.names(gp@tab),]

Dgeo = distm(x = site_dat.ordered[,c("lon", "lat")], fun = distVincentyEllipsoid)
rownames(Dgeo) = site_dat.ordered$state
colnames(Dgeo) = site_dat.ordered$state
Dgeo = Dgeo/1000

mantel(Dgen.2, Dgeo)
# Mantel statistic r: 0.1575 
# Significance: 0.203 

# about the same as the above