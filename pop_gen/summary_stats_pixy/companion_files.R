library(dplyr)
library(vcfR)

vcf <- read.vcfR("data/Nf/final_tables/rm_dups/FINAL_snp.IBD_analyses.vcf.gz", verbose = FALSE)
sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")

#pops file is tab-separated sample_name popID

sample_metadata.Nf %>% head

write.table((sample_metadata.Nf %>% select(Sequence_label, state)), "data/sample_metadata/Nf.pixy_pops.tsv", sep = "\t", col.names = F, row.names = F, quote = F)


state_n = sample_metadata.Nf %>% group_by(state) %>% summarise(n = n())
low_n = state_n %>% filter(n < 4)
low_n

write.table((sample_metadata.Nf %>% filter(!state %in% low_n$state) %>% select(Sequence_label, state)), "data/sample_metadata/Nf.tajD_pops.tsv", sep = "\t", col.names = F, row.names = F, quote = F)


            