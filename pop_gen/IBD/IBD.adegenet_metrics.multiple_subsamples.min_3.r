require(vcfR)
require(tidyr)
require(dplyr)
require(ggplot2)
require(pegas)
require(reshape2)
require(adegenet)
#require(dartR) #no longer needed bc not using gl2gi
source("R_scripts/ggplot_theme.txt")
set.seed(12345)

#metadata
source("R_scripts/make_site_metadata.r")
site_coords = read.table("data/sample_metadata/site_coords_for_map.txt", header = T)
ind.metrics = left_join(data.frame(sample = fam_info[,1]), sample_metadata %>% select(sample,state.name) ) %>%
left_join(., site_coords %>% select(state.name, lat, lon))
#########


#filtered VCF
vcf <- read.vcfR("data/Nf_SPANDx_all_seqs/out.filtered.LD_filtered_0.5_10Kb.vcf", verbose = FALSE)
gl = vcfR2genlight(vcf)
#filter to the remaining samples
ind.metrics.subset = ind.metrics[ind.metrics$sample %in% gl@ind.names,]

#Set metadata in genLight
gl@other$ind.metrics = ind.metrics.subset
gl@other$latlong = ind.metrics.subset[,3:4]
gl@other$latlong
gl@pop = as.factor(ind.metrics.subset$state.name) #need to set pop for the ibd test to work
gl@ploidy = rep(as.integer(1), nInd(gl))
nPop(gl)
pop(gl)

#SOME SITES HAVE SMALL SAMPLE SIZE
#Set min sample size to 3
low_n_states = ( (gl@other$ind.metrics %>% group_by(state.name) %>% summarize(n = n()) ) %>% filter(n < 3) )$state.name
keep.ind.list = data.frame(gl@ind.names, pop.gl = pop(gl)) %>% filter(!pop.gl %in% c(low_n_states))

#loop through this routine twenty times (i.e., twenty random samples) and calculate dissimilarity matrices, THEN average matrices at the end.
#will keep dissim matrices in a list

distances.list = list()
#based on the lowest number of samples to include sites
min_samps = 3

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
saveRDS(distances.list, "data/intermediate_RDS/DSCE.three_samples_per_site.rds")
distances.list = readRDS("data/intermediate_RDS/DSCE.three_samples_per_site.rds")
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

ibd.2 <- mantel.randtest(mean_Dgen,Dgeo) # r = -0.04825996  , P = 0.569
ibd.2.log <- mantel.randtest(mean_Dgen,log(Dgeo)) # r = -0.006581043 , P = 0.463
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

pdf("figures/IBD.Dcse_test.3_samples.pdf", width = 9, height = 5)
p1
p2
dev.off()

######################################
#Running without VA (very young site)#
######################################


mean_Dgen.no_VA = as.matrix(mean_Dgen)[-13,-13] %>% as.dist()
Dgeo.no_VA = as.matrix(Dgeo)[-13,-13] %>% as.dist()

#test
ibd.2 <- mantel.randtest(mean_Dgen.no_VA,Dgeo.no_VA) # r = -0.1270664   , P = 0.414
ibd.2.log <- mantel.randtest(mean_Dgen.no_VA,log(Dgeo.no_VA)) # r = -0.1043556 , P = 0.697
print(ibd.2)
print(ibd.2.log)

#note that the log is better justified based on distance decay and 2-D spread (refs)

