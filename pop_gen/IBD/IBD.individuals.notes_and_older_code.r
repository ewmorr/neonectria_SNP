library(vcfR)
#library(tidyr)
#library(dplyr)
library(ggplot2)
#library(pegas)
#library(reshape2)
library(adegenet)

source("library/ggplot_theme.txt")
#set.seed(12345)

#metadata
sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
sample_metadata.Nd = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
#########

#########
#########
#We were originally using adegenet to calculate populations level distances
# We would like to instead compute pairwise distances between inds
# The genind format can not handle poly alleles (only bialleles)
# The most basic method to calculate distance is using dist() euclidean
# but this fails for polyalleles because different character states will be 
# treated as greater than simple REF:ALT (e.g., 1-0 versus 2-0 vs 2-1).
# Can calculate pairwise Hamming distance (character difs) using 
# phangorn::dist.dna with pairwise deletion and using the nt fasta written
# from the VCF (as used for the phylogeny). Use model = "raw" to compute ratio
# and then multiply by the number of sites for count based Hamming (which is
# suitable for poppr::poppr.amova, e.g.)
# 
# phangorn::dist.dna(x = fasta, model = "raw", pairwise.deletion = T)
# OR
# phangorn::dist.hamming(x = fasta, ratio = T, exclude = "pairwise")


#filtered VCF
vcf.Nf <- read.vcfR("data/Nf/final_tables/rm_dups/FINAL_snp.IBD_analyses.vcf.gz", verbose = FALSE)
gl.Nf = vcfR2genlight(vcf.Nf) #note that genlight object only supports bialleles
# 47047 alleles excluded
# #note that this does not match the number of sites removed by vcftools (666201)
rm(vcf.Nf)
gc()

vcf.Nd <- read.vcfR("data/Nd/final_tables/rm_dups/FINAL_snp.IBD_analyses.vcf.gz", verbose = FALSE)
gl.Nd = vcfR2genlight(vcf.Nd)
rm(vcf.Nd)
gc()


#Set metadata in genLight
gl@other$ind.metrics = sample_metadata.Nf



gl@other$latlong = ind.metrics.subset[,3:4]
gl@other$latlong
gl@pop = as.factor(ind.metrics.subset$state.name) #need to set pop for the ibd test to work
gl@ploidy = rep(as.integer(1), nInd(gl))
nPop(gl)
pop(gl)

#SOME SITES HAVE SMALL SAMPLE SIZE
#Set min sample size to 5
low_n_states = ( (gl@other$ind.metrics %>% group_by(state.name) %>% summarize(n = n()) ) %>% filter(n < 5) )$state.name
keep.ind.list = data.frame(gl@ind.names, pop.gl = pop(gl)) %>% filter(!pop.gl %in% c(low_n_states))

#loop through this routine twenty times (i.e., twenty random samples) and calculate dissimilarity matrices, THEN average matrices at the end.
#will keep dissim matrices in a list

distances.list = list()
#based on the lowest number of samples to include sites
min_samps = 5

for(u in 1:100){
    
    keep.ind.rand = list()
    pops_incl = keep.ind.list$pop.gl %>% as.character %>% unique
    
    #select min_samps random samples from each site for distance calc
    for(i in 1:length(pops_incl)){
        temp = keep.ind.list %>% filter(pop.gl == pops_incl[i])
        keep.ind.rand[[pops_incl[i]]] = temp[sample(1:nrow(temp), min_samps), ] # random sample of min_samps rows
    }

    keep.ind.rand.df = bind_rows(keep.ind.rand)
    ##########################################

    gl.subset = gl[gl@ind.names %in% keep.ind.rand.df$gl.ind.names]
    
    #Need to convert to genpop for adegenet dist.genpop
    #INSTEAD OF DARTR METHOD TRY THE DATAFRAME CONVERSION AS RECOMMENDED BY ADEGENET AUTHORS (https://lists.r-forge.r-project.org/pipermail/adegenet-forum/2014-May/000840.html)
    #First convert to a data.frame whih will give a table of 0, 1, NA, and then add 1 to values to have correct conversion of NA (0 is default)
    y = as.data.frame(gl.subset)
    y = y + 1
    #y[is.na(y)] = 0
    #colnames(y) = 1:ncol(y) #this is needed when using the fasta conversion instead of vcf input bc the colnames get funny
    #colnames(y)
    gi.subset = df2genind(y, ploidy=1)
    #reset pop
    gi.subset@pop = gl.subset@pop

    #Remove uninformative sites
    #retrieve the colnames of sites with only one allele
    to_remove = names(gi.subset@loc.n.all[gi.subset@loc.n.all == 1 ])
    if(length(to_remove) > 0){
        #get the col index
        rm_indx = which(colnames(y) %in% to_remove)
        gi.subset.rm = gi.subset[loc=-rm_indx]

        gp = genind2genpop(gi.subset.rm)
    }else{
        gp = genind2genpop(gi.subset)
    }

    #There are several distance metrics available
    #Method 2 is "Angular distance or Edward's distance" D[CSE]
    Dgen.2 <- dist.genpop(gp,method=2)
    distances.list[[u]] = Dgen.2
}

saveRDS(distances.list, "data/intermediate_RDS/DSCE.five_samples_per_site.rds")
distances.list = readRDS("data/intermediate_RDS/DSCE.five_samples_per_site.rds")

#average gen dist matrix
mean_Dgen = Reduce("+", distances.list) / length(distances.list)

#Get pariwise site distances

#filter to included sites
site_coords.subset = site_coords %>% filter(state.name %in% levels(pop(gl.subset)))
#need to order the df to compare distances correctly
site_order = levels(pop(gl.subset))
site_coords.subset.order = site_coords.subset[match(site_order, site_coords.subset$state.name), ]
rownames(site_coords.subset.order) = site_coords.subset.order$state.name
#calculate Mercator prjected distance
Dgeo <- dist(dismo::Mercator(site_coords.subset.order[,c("lon", "lat")]))
saveRDS(Dgeo, "data/intermediate_RDS/geo_distance.min_5_samples.rds")

ibd.2 <- mantel.randtest(mean_Dgen,Dgeo) # r = 0.2191161  , P = 0.137
ibd.2.log <- mantel.randtest(mean_Dgen,log(Dgeo)) # r = 0.2700557 , P = 0.099
print(ibd.2)
print(ibd.2.log)

#######################
#Reformat for plotting#
#for plotting
#######################
Dgeo.long = subset(reshape2::melt(Dgeo %>% as.matrix), value != 0)
Dgen.long = subset(reshape2::melt(mean_Dgen %>% as.matrix), value != 0)

Dgen.Dgeo = data.frame(gen = Dgen.long$value, geo = Dgeo.long$value)

p1 = ggplot(Dgen.Dgeo, aes(x = geo/1000, y = gen)) +
geom_point() +
geom_smooth(method = "lm", color = "black", se = F, linetype = 2) +
labs(x = "Pairwise site geographic distance (km)", y = expression(paste("Pairwise site genetic distance (D"["CSE"], ")"))) +
my_gg_theme

p2 = ggplot(Dgen.Dgeo, aes(x = geo/1000, y = gen)) +
geom_point() +
geom_smooth(method = "lm", color = "black", se = F, linetype = 2) +
labs(x = "Pairwise site geographic distance (km)", y = expression(paste("Pairwise site genetic distance (D"["CSE"], ")"))) +
scale_x_continuous(trans = "log", breaks = c(200, 500, 1000, 2000)) +
my_gg_theme

pdf("figures/IBD.Dcse_test.pdf", width = 9, height = 5)
p1
p2
dev.off()

######################################
#Running without VA (very young site)#
######################################


mean_Dgen.no_VA = as.matrix(mean_Dgen)[-11,-11] %>% as.dist()
Dgeo.no_VA = as.matrix(Dgeo)[-11,-11] %>% as.dist()

#test
ibd.2 <- mantel.randtest(mean_Dgen.no_VA,Dgeo.no_VA) # r = 0.4338682  , P = 0.06
ibd.2.log <- mantel.randtest(mean_Dgen.no_VA,log(Dgeo.no_VA)) # r = 0.4936303 , P = 0.012
print(ibd.2)
print(ibd.2.log)

#note that the log is better justified based on distance decay and 2-D spread (refs)

#stats from the previous run before repolishing the genome and including a few extra sample
#no log
# r =  0.0.4907064 , P = 0.027
#log
# r = 0.5511402  , P = 0.01

#######################
#Reformat for plotting#
#for plotting
#######################
Dgeo.long = subset(reshape2::melt(Dgeo.no_VA %>% as.matrix), value != 0)
Dgen.long = subset(reshape2::melt(mean_Dgen.no_VA %>% as.matrix), value != 0)

Dgen.Dgeo = data.frame(gen = Dgen.long$value, geo = Dgeo.long$value)

p1 = ggplot(Dgen.Dgeo, aes(x = geo/1000, y = gen)) +
geom_point() +
geom_smooth(method = "lm", color = "black", se = F, linetype = 2) +
labs(x = "Pairwise site geographic distance (km)", y = expression(paste("Pairwise site genetic distance (D"["CSE"], ")"))) +
my_gg_theme

p2 = ggplot(Dgen.Dgeo, aes(x = geo/1000, y = gen)) +
geom_point() +
geom_smooth(method = "lm", color = "black", se = F, linetype = 2) +
labs(x = "Pairwise site geographic distance (km)", y = expression(paste("Pairwise site genetic distance (D"["CSE"], ")"))) +
scale_x_continuous(trans = "log", breaks = c(200, 500, 1000, 2000)) +
my_gg_theme


pdf("figures/IBD.Dcse_test.no_VA.pdf", width = 9, height = 5)
p1
p2
dev.off()
