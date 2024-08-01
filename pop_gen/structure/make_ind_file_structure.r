library(vcfR)

#just pulling in dp to get sample IDs for ordering
vcf <- read.vcfR("data/Nf/final_tables/rm_dups/FINAL_snp.mac_ge2.biallele.LD.structure_analyses.vcf.gz", verbose = FALSE)
sample_ids <- extract.gt(vcf, element='DP', as.numeric=TRUE) %>% colnames()

metadata = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")

ind_file = left_join(data.frame(Sequence_label = sample_ids), metadata %>% dplyr::select(Sequence_label, state) )
write.table(ind_file, "pop_gen/structure/Nf_ind_file.structure", quote = F, row.names = F, col.names = F)

#Nd
vcf <- read.vcfR("data/Nd/final_tables/rm_dups/FINAL_snp.mac_ge2.biallele.LD.structure_analyses.vcf.gz", verbose = FALSE)
sample_ids <- extract.gt(vcf, element='DP', as.numeric=TRUE) %>% colnames()

metadata = read.csv("data/sample_metadata/Nd_filtered.lat_lon_dur_inf.csv")

ind_file = left_join(data.frame(Sequence_label = sample_ids), metadata %>% dplyr::select(Sequence_label, state) )
write.table(ind_file, "pop_gen/structure/Nd_ind_file.structure", quote = F, row.names = F, col.names = F)
