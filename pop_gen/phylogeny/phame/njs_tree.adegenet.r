require(adegenet)
require(ape)
require(phangorn)
library(dplyr)


#read tree
tre = read.tree(file = "data/phame_phylogeny_no_trim/results/trees/core_fugr_all_snp_alignment.fna.treefile")
tre
str(tre)
tre$tip.label
tre$tip.label = sub("_pread", "", tre$tip.label)
#the ref genomes are indices 145, 146, 147

plot(tre)
plot(tre, type = "unrooted", show.tip = F)
tre = drop.tip(tre, match(c("nedi", "neco", "nefa"), tre$tip.label))
plot(tre)
#NG149, NG20, NG152, NG27 need to be removed as well as these are the dups
tre = drop.tip(tre, match(c("NG149", "NG20", "NG152", "NG27"), tre$tip.label))
plot(tre)



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

#root on fugr
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