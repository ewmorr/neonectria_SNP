set.seed(12345)
library(vcfR)
library(dplyr)
library(adegenet)

#premise conda env
#mamba create --name R-pop_gen conda-forge::r-tidyr conda-forge::r-dplyr bioconda::r-adegenet bioconda::r-vcfr
#conda activate R-pop_gen

#metadata
sample_metadata.Nf = read.csv("~/Nf_pop_IBD_11182024/Nf_filtered.lat_lon_dur_inf.csv")


#########
#prev used table of sample,state.name,lat,lon
ind.metrics = sample_metadata.Nf %>% select(Sequence_label, state, lat, lon)

#filtered VCF
vcf <- read.vcfR("~/Nf_pop_IBD_11182024/FINAL_snp.IBD_analyses.vcf.gz", verbose = FALSE)
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

#we previously picked 3 samples at random without replacement
#instead we will bootstrap the mean using n of three to account for unequal 
# sample size
# then average matrices at the end.
#will keep dissim matrices in a list

distances.list = list()
#based on the lowest number of samples to include sites
min_samps = 3

n_boots = 1000
#run the bootstrapping
for(u in 1:n_boots){
    
    keep.ind.rand = list()
    pops_incl = keep.ind.list$pop.gl %>% as.character %>% unique
    
    #select min_samps random samples from each site for distance calc
    for(i in 1:length(pops_incl)){
        temp = keep.ind.list %>% filter(pop.gl == pops_incl[i])
        keep.ind.rand[[pops_incl[i]]] = temp[sample(1:nrow(temp), min_samps, replace = T),] # random sample of min_samps rows
    }

    keep.ind.rand.df = bind_rows(keep.ind.rand)
    ##########################################

    gl.subset = gl[gl@ind.names %in% keep.ind.rand.df$gl.ind.names]
    
    #Need to convert to genpop for adegenet dist.genpop
    
    #First convert to a data.frame which will give a table of 0, 1, NA, and then add 1 to values to have correct conversion of NA (0 is default)
    y = as.data.frame(gl.subset)
    gi.subset = df2genind(y + 1, ploidy=1)
    #reset pop
    gi.subset@pop = gl.subset@pop

    #Remove uninformative sites
    #retrieve the colnames of sites with only one allele
    to_remove = names(gi.subset@loc.n.all[gi.subset@loc.n.all == 1 ])
    #get the col index
    rm_indx = which(colnames(y) %in% to_remove)
    gi.subset.rm = gi.subset[loc=-rm_indx]

    gp = genind2genpop(gi.subset.rm)

    #There are several distance metrics available
    #Method 2 is "Angular distance or Edward's distance" D[CSE]
    Dgen.2 <- dist.genpop(gp,method=2)
    distances.list[[u]] = Dgen.2
}

saveRDS(distances.list, "~/Nf_pop_IBD_11182024/Nf.DSCE.three_samples_per_site.rds")
