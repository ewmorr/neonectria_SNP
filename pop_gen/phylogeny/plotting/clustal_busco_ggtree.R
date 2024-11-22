library(dplyr)
library(ggplot2)
library(ggtree)

tre = treeio::read.tree("data/shared_buscos/final_tables/rm_dups/FINAL_snp.snps_only.for_phylogeny.ph")
str(tre)


ggtree(tre)
