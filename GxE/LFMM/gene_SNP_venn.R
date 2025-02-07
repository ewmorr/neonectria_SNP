
library(ggVennDiagram)
library(ggplot2)
library(dplyr)

#read gene IDs
ppt_genes = read.table("data/Nf_LFMM_tables/ppt.geneIDs.nearest_neighbors.txt")
hdd4_genes = read.table("data/Nf_LFMM_tables/hdd4.geneIDs.nearest_neighbors.txt")
freezeThaw_genes = read.table("data/Nf_LFMM_tables/freezeThaw.geneIDs.nearest_neighbors.txt")
tmin_genes = read.table("data/Nf_LFMM_tables/tmin.geneIDs.nearest_neighbors.txt")
growing_hdd4_genes = read.table("data/Nf_LFMM_tables/growing_hdd4.geneIDs.nearest_neighbors.txt")

#5 var Venn
gene_list = list(
    "nongrowing\nGDD4" = hdd4_genes$V1,
    "growing\nGDD4" = growing_hdd4_genes$V1,
    "freeze-thaw" =  freezeThaw_genes$V1,
    "Tmin" =  tmin_genes$V1,
    "MAP" = ppt_genes$V1
)

ggVennDiagram(gene_list) + 
    scale_fill_gradient(low="blue",high = "red") +
    scale_color_manual(values = rep("grey", 5) )

p1 = ggVennDiagram(gene_list) + 
    scale_fill_gradient(low="blue",high = "red") +
    scale_color_manual(values = rep("grey", 5) )

pdf("figures/LFMM_genes_venn.5_vars.pdf")
p1
dev.off()

#################################
#nongrowing_GDD, MAP, freeze-thaw
gene_list = list(
    "nongrowing\nGDD4" = hdd4_genes$V1,
    "freeze-thaw" =  freezeThaw_genes$V1,
    "MAP" = ppt_genes$V1
)

ggVennDiagram(gene_list) + 
    scale_fill_gradient(low="blue",high = "red") +
    scale_color_manual(values = rep("grey", 5) )

p1 = ggVennDiagram(gene_list) + 
    scale_fill_gradient(low="blue",high = "red") +
    scale_color_manual(values = rep("grey", 5) )

pdf("figures/LFMM_genes_venn.nongrowing_hdd4.tmin.freezeThaw.pdf", width = 4.5, height = 4.5)
p1
dev.off()

#######################
#Output lists of gene IDs shared between vars and unique to vars
#Using the three var list
#Tmin and growing season GDD4 have only 1-2 genes each
#and both end up in the intersection of nongrowing GDD4 qnd MAP

#conveniently all ft genes are unique
write.table(freezeThaw_genes, "data/Nf_LFMM_tables/freezeThaw.geneIDs.unique.txt", col.names = F, row.names = F, quote = F)

hdd4_ppt_shared = ppt_genes %>% filter(V1 %in% hdd4_genes$V1)
ppt_unique = ppt_genes %>% filter(!V1 %in% hdd4_genes$V1)
hdd4_unique = hdd4_genes %>% filter(!V1 %in% ppt_genes$V1)

write.table(hdd4_ppt_shared, "data/Nf_LFMM_tables/hdd4_ppt_shared.geneIDs.shared.txt", col.names = F, row.names = F, quote = F)
write.table(ppt_unique, "data/Nf_LFMM_tables/ppt.geneIDs.unique.txt", col.names = F, row.names = F, quote = F)
write.table(hdd4_unique, "data/Nf_LFMM_tables/hdd4.geneIDs.unique.txt", col.names = F, row.names = F, quote = F)

################
#temp vars
gene_list = list(
    "nongrowing\nGDD4" = hdd4_genes$V1,
    "growing\nGDD4" = growing_hdd4_genes$V1,
    "freeze-thaw" =  freezeThaw_genes$V1,
    "Tmin" =  tmin_genes$V1
)

ggVennDiagram(gene_list) + 
    scale_fill_gradient(low="blue",high = "red") +
    scale_color_manual(values = rep("grey", 5) )

p1 = ggVennDiagram(gene_list) + 
    scale_fill_gradient(low="blue",high = "red") +
    scale_color_manual(values = rep("grey", 5) )

pdf("figures/LFMM_genes_venn.temp_vars.pdf")
p1
dev.off()

###############################################
#SNP Venn diagrams
##################

ppt_snps = read.table("data/Nf_LFMM_tables/ppt_lfmm.txt", header = T)
hdd4_snps = read.table("data/Nf_LFMM_tables/hdd4_lfmm.txt", header = T)
growing_hdd4_snps = read.table("data/Nf_LFMM_tables/growing_hdd4_lfmm.txt", header = T)
tmin_snps = read.table("data/Nf_LFMM_tables/tmin_lfmm.txt", header = T)
freezeThaw_snps = read.table("data/Nf_LFMM_tables/freezeThaw_lfmm.txt", header = T)

#filter to significant and paste name column
ppt_snps.sig = ppt_snps %>% filter(FDR.sig == "sig") %>% select(scaffold, position)
ppt_snps.sig$name = paste(ppt_snps.sig$scaffold, ppt_snps.sig$position, sep = "-")

hdd4_snps.sig = hdd4_snps %>% filter(FDR.sig == "sig") %>% select(scaffold, position)
hdd4_snps.sig$name = paste(hdd4_snps.sig$scaffold, hdd4_snps.sig$position, sep = "-")

growing_hdd4_snps.sig = growing_hdd4_snps %>% filter(FDR.sig == "sig") %>% select(scaffold, position)
growing_hdd4_snps.sig$name = paste(growing_hdd4_snps.sig$scaffold, growing_hdd4_snps.sig$position, sep = "-")

tmin_snps.sig = tmin_snps %>% filter(FDR.sig == "sig") %>% select(scaffold, position)
tmin_snps.sig$name = paste(tmin_snps.sig$scaffold, tmin_snps.sig$position, sep = "-")

freezeThaw_snps.sig = freezeThaw_snps %>% filter(FDR.sig == "sig") %>% select(scaffold, position)
freezeThaw_snps.sig$name = paste(freezeThaw_snps.sig$scaffold, freezeThaw_snps.sig$position, sep = "-")



#5 var Venn
snp_list = list(
    "nongrowing\nGDD4" = hdd4_snps.sig$name,
    "growing\nGDD4" = growing_hdd4_snps.sig$name,
    "freeze-thaw" =  freezeThaw_snps.sig$name,
    "Tmin" =  tmin_snps.sig$name,
    "MAP" = ppt_snps.sig$name
)

ggVennDiagram(snp_list) + 
    scale_fill_gradient(low="blue",high = "red") +
    scale_color_manual(values = rep("grey", 5) )

p1 = ggVennDiagram(snp_list) + 
    scale_fill_gradient(low="blue",high = "red") +
    scale_color_manual(values = rep("grey", 5) )

pdf("figures/LFMM_snps_venn.5_vars.pdf")
p1
dev.off()

#nongrowing_hdd4.tmin.freezeThaw
snp_list = list(
    "nongrowing\nGDD4" = hdd4_snps.sig$name,
    "freeze-thaw" =  freezeThaw_snps.sig$name,
    "MAP" = ppt_snps.sig$name
)

ggVennDiagram(snp_list) + 
    scale_fill_gradient(low="blue",high = "red") +
    scale_color_manual(values = rep("grey", 5) )

p1 = ggVennDiagram(snp_list) + 
    scale_fill_gradient(low="blue",high = "red") +
    scale_color_manual(values = rep("grey", 5) )

pdf("figures/LFMM_snps_venn.nongrowing_hdd4.tmin.freezeThaw.pdf", width = 4.5, height = 4.5)
p1
dev.off()

#temp vars
snp_list = list(
    "nongrowing\nGDD4" = hdd4_snps.sig$name,
    "growing\nGDD4" = growing_hdd4_snps.sig$name,
    "freeze-thaw" =  freezeThaw_snps.sig$name,
    "Tmin" =  tmin_snps.sig$name
)

ggVennDiagram(snp_list) + 
    scale_fill_gradient(low="blue",high = "red") +
    scale_color_manual(values = rep("grey", 5) )

p1 = ggVennDiagram(snp_list) + 
    scale_fill_gradient(low="blue",high = "red") +
    scale_color_manual(values = rep("grey", 5) )

pdf("figures/LFMM_snps_venn.temp_vars.pdf")
p1
dev.off()
