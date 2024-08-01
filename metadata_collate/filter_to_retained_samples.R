library(vcfR)
library(dplyr)

# NF

#filtered VCF
vcf <- read.vcfR("data/Nf/final_tables/rm_dups/FINAL_snp.mac_ge2.biallele.LD.structure_analyses.vcf.gz", verbose = FALSE)

#just pulling in dp to get sample IDs
sample_ids <- extract.gt(vcf, element='DP', as.numeric=TRUE) %>% colnames()

rm(vcf)
gc()

metadata = read.csv("data/sample_metadata/lat_lon_dur_inf.csv")

metadata.Nf = metadata[metadata$Sequence_label %in% sample_ids,]
nrow(metadata.Nf)

#there are some samples of different individuals from the same tree; 
# N149, 118 (ANF1 1); NG114, 152 (ANF1 10); 
#look at which ones have better completeness
# run with vcftools on original table in CLI
# 
#NG118   1116217 0       24087   0.0215791
#NG149   1116217 0       45068   0.0403757
#NG152   1116217 0       25970   0.0232661
#NG114   1116217 0       10974   0.00983142
metadata.Nf = metadata.Nf[!metadata.Nf$Sequence_label %in% c("NG149", "NG152"),]
nrow(metadata.Nf)

write.csv(metadata.Nf, "data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv", row.names = F)

# Nd

#filtered VCF
vcf <- read.vcfR("data/Nd/final_tables/FINAL_snp.mac_ge2.biallele.LD.strctr_analyses.vcf.gz", verbose = FALSE)

#just pulling in dp to get sample IDs
sample_ids <- extract.gt(vcf, element='DP', as.numeric=TRUE) %>% colnames()
rm(vcf)
gc()

metadata.Nd = metadata[metadata$Sequence_label %in% sample_ids,]
nrow(metadata.Nd)

#there are some samples of different individuals from the same tree or the same indv sequenced twice; 
# NG163, NG20 (MI1 1.1.1); NG27, NG144 (MI1 2); 
#look at which ones have better completeness
# run with vcftools on original table in CLI
# 
#NG20    1600316 0       62002   0.0387436
#NG27    1600316 0       46514   0.0290655
#NG144   1600316 0       18411   0.0115046
#NG163   1600316 0       13410   0.0083796
metadata.Nd = metadata.Nd[!metadata.Nd$Sequence_label %in% c("NG20", "NG27"),]
nrow(metadata.Nd)

write.csv(metadata.Nd, "data/sample_metadata/Nd_filtered.lat_lon_dur_inf.csv", row.names = F)
