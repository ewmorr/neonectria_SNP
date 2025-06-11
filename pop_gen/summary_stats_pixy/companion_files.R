library(dplyr)
library(vcfR)

vcf <- read.vcfR("data/Nf/final_tables/rm_dups/FINAL_snp.IBD_analyses.vcf.gz", verbose = FALSE)
sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")

#pops file is tab-separated sample_name popID

sample_metadata.Nf %>% head

write.table((sample_metadata.Nf %>% select(Sequence_label, state)), "data/sample_metadata/Nf.pixy_pops.tsv", sep = "\t", col.names = F, row.names = F, quote = F)

write.table( data.frame(label = sample_metadata.Nf %>% select(Sequence_label), pop = rep("a", nrow(sample_metadata.Nf)) ), "data/sample_metadata/Nf.single_pixy_pop.tsv", sep = "\t", col.names = F, row.names = F, quote = F)

state_n = sample_metadata.Nf %>% group_by(state) %>% summarise(n = n())
low_n = state_n %>% filter(n < 2)
low_n

write.table((sample_metadata.Nf %>% filter(!state %in% low_n$state) %>% select(Sequence_label, state)), "data/sample_metadata/Nf.tajD_pops.tsv", sep = "\t", col.names = F, row.names = F, quote = F)

pca_clust = read.csv("data/sample_metadata/Nf_filtered.HCPC_clust.csv")
pca_clust

pca_clust$HCPC.cluster = paste0("clust_", pca_clust$HCPC.cluster)
    
write.table(pca_clust, "data/sample_metadata/Nf_pixy_pops.HCPC_clust.txt", sep = "\t", col.names = F, row.names = F, quote = F)
    

# Nd
vcf <- read.vcfR("data/Nd/final_tables/rm_dups/FINAL_snp.IBD_analyses.vcf.gz", verbose = FALSE)
sample_metadata.Nd = read.csv("data/sample_metadata/Nd_filtered.lat_lon_dur_inf.csv")

write.table( data.frame(label = sample_metadata.Nd %>% select(Sequence_label), pop = rep("a", nrow(sample_metadata.Nd)) ), "data/sample_metadata/Nd.single_pixy_pop.tsv", sep = "\t", col.names = F, row.names = F, quote = F)

# Nc
vcf <- read.vcfR("data/Nc/final_tables/rm_dups/FINAL_snp.IBD_analyses.vcf.gz", verbose = FALSE)
sample_metadata.Nc = read.csv("data/sample_metadata/Nc_canton_loc_date.lat_lon.csv")

write.table( data.frame(label = sample_metadata.Nc %>% select(Sequence_label), pop = rep("a", nrow(sample_metadata.Nc)) ), "data/sample_metadata/Nc.single_pixy_pop.tsv", sep = "\t", col.names = F, row.names = F, quote = F)
