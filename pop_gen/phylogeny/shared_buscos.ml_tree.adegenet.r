require(adegenet)
require(ape)
require(phangorn)



# read SNP multi seq aligment fasta
seqs = fasta2DNAbin(file = "data/shared_buscos/final_tables/rm_dups/FINAL_snp.snps_only.for_phylogeny.fasta")
#695,811 nt

# https://adegenet.r-forge.r-project.org/files/PRstats/practical-introphylo.1.0.pdf
# convert to phangorn format
seqs.phyDat = as.phyDat(seqs)

#initialize tree with NJ (using original DNAbinobject, not phangorn format)
#
dist_mat = dist.dna(seqs, pairwise.deletion = T, model="raw") #must use pairwise deletion because of high number of Ns
sum(is.na(dist_mat) )
sum(!is.na(dist_mat) )
range_dist = range(dist_mat, na.rm = T)
range_dist[2]/range_dist[1]
# the TN93 model is giving NAs. The docs say that most evolutionary distances are 
# undefined with a lot of difs (e.g., above 0.75). We return values up to 0.825
# with the raw proportion.
# 

tre.ini = nj(dist_mat) 
plot(tre.ini)

#optimization
fit.ini = pml(tre.ini, seqs.phyDat, k = 4)

#this takes quite a while to run (like 30 minutes or more)
fit = optim.pml(fit.ini, optNni=T, optBf = T, optQ = T, optGamma = T)

# computed
fit

fit$tree
plot(fit$tree)
anova(fit.ini,fit)

#bootstrap.pml() for ml tree or boot.phylo() for NJ
## this paper has workflows https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11117635/

#write the ml tree before plotting so don't have to rerun

saveRDS(fit, "data/shared_buscos/phylogeny/ML_tree.rds")
