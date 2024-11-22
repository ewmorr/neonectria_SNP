require(adegenet)
require(ape)
require(phangorn)
library(dplyr)



# read SNP multi seq aligment fasta
seqs = fasta2DNAbin(file = "data/Fugr1_ref/core.Fusgr1-neonectria.snps_plus_invariant_aln.fasta")
#290,590 nt

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
tre = bionj(dist_mat) 
plot(tre)
plot(tre, type = "unrooted", show.tip = F)

#Read metadata
#
sample_metadata = read.csv("data/sample_metadata/core_fugr.lat_lon_dur_inf.csv")
colnames(sample_metadata)
colnames(sample_metadata)[1] = "label" #for join below
annot = left_join(data.frame(label = tre$tip.label), sample_metadata)
annot$color = vector(mode = "character", length = nrow(annot))
annot[annot$spp == "Nf", "color"] = "blue"
annot[annot$spp == "Nc", "color"] = "red"
annot[annot$spp == "Nd", "color"] = "green"
annot[annot$spp == "Fg", "color"] = "yellow"

####################################
plot(tre, type = "unrooted", show.tip = F)
tiplabels(annot$spp, bg=annot$color,
          cex=.5)

#root on Neco
tre2 = root(tre, out = 1)
tre2 = ladderize(tre2)
plot(tre2, , show.tip = F)
tiplabels(annot$spp, bg=annot$color,
          cex=.5)

#root on Nedi
tre2 = root(tre, out = 103) 
tre2 = ladderize(tre2)
plot(tre2, , show.tip = F)
tiplabels(annot$spp, bg=annot$color,
          cex=.5)



############################
#############################