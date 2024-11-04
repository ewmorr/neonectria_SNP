library(vcfR)
library(dplyr)

#thinned Nf
vcf <- read.vcfR("data/Nf/final_tables/rm_dups/FINAL_snp.mac_ge2.biallele.LD.structure_analyses.vcf.gz", verbose = FALSE)
gt <- extract.gt(vcf, element='GT', as.numeric=TRUE)
head(gt)
nrow(gt)
NA_prop = rowSums(is.na(gt))/ncol(gt)
plot(sort(NA_prop))

NA_prop_df = data.frame(
    index = 1:length(NA_prop),
    NA_prop = NA_prop
)
NA_prop_df$col = ifelse(NA_prop_df$NA_prop == 0, "red", "black")

no_NA_names = names(NA_prop[NA_prop == 0])
length(no_NA_names)
#20465
gt.no_na = gt[row.names(gt) %in% no_NA_names,]
nrow(gt.no_na)
gt = gt.no_na
gt[is.na(gt)] = -9
head(gt)
gt_str = t(gt)
nrow(gt_str)
#115 individuals (for structure input options)
n_snps = ncol(gt_str)
n_snps
#20465 snps  (for structure input options)
colnames(gt_str) = paste("SNP_", seq(1, n_snps), sep = "")
head(gt_str[,1:5])

#adding the loc info
sample_ids = row.names(gt_str)
metadata = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
loc_data = left_join(data.frame(Sequence_label = sample_ids), metadata %>% dplyr::select(Sequence_label, state) )
loc_int = data.frame(state = unique(loc_data$state), int = 1:length(unique(loc_data$state)))
loc_data.int = left_join(loc_data, loc_int, by = "state")
gt_str.loc = data.frame(LocData = loc_data.int$int, gt_str)
gt_str.loc[,1:5] %>% head

write.table(gt_str.loc, "data/Nf/final_tables/rm_dups/FINAL_snp.thinned_no_NA.LocData.structure", quote = F, sep = "\t", row.names = T, col.names = NA)
#subsequently delete the LocData colname manually or with sed

# Nd
vcf <- read.vcfR("data/Nd/final_tables/rm_dups/FINAL_snp.mac_ge2.biallele.LD.structure_analyses.vcf.gz", verbose = FALSE)
gt <- extract.gt(vcf, element='GT', as.numeric=TRUE)
head(gt)
nrow(gt)
NA_prop = rowSums(is.na(gt))/ncol(gt)
plot(sort(NA_prop))

NA_prop_df = data.frame(
    index = 1:length(NA_prop),
    NA_prop = NA_prop
)
NA_prop_df$col = ifelse(NA_prop_df$NA_prop == 0, "red", "black")

no_NA_names = names(NA_prop[NA_prop == 0])
length(no_NA_names)
#57579
gt.no_na = gt[row.names(gt) %in% no_NA_names,]
nrow(gt.no_na)
gt = gt.no_na
gt[is.na(gt)] = -9
head(gt)
gt_str = t(gt)
nrow(gt_str)
#30 individuals (for structure input options)
n_snps = ncol(gt_str)
n_snps
#57579 snps  (for structure input options)
colnames(gt_str) = paste("SNP_", seq(1, n_snps), sep = "")
head(gt_str[,1:5])

#adding the loc info
sample_ids = row.names(gt_str)
metadata = read.csv("data/sample_metadata/Nd_filtered.lat_lon_dur_inf.csv")
loc_data = left_join(data.frame(Sequence_label = sample_ids), metadata %>% dplyr::select(Sequence_label, state) )
loc_int = data.frame(state = unique(loc_data$state), int = 1:length(unique(loc_data$state)))
loc_data.int = left_join(loc_data, loc_int, by = "state")
gt_str.loc = data.frame(LocData = loc_data.int$int, gt_str)
gt_str.loc[,1:5] %>% head

write.table(gt_str.loc, "data/Nd/final_tables/rm_dups/FINAL_snp.thinned_no_NA.LocData.structure", quote = F, sep = "\t", row.names = T, col.names = NA)
#subsequently delete the LocData colname manually or with sed



# Nc
vcf <- read.vcfR("data/Nc/final_tables/FINAL_snp.mac_ge2.biallele.LD.structure_analyses.vcf.gz", verbose = FALSE)
gt <- extract.gt(vcf, element='GT', as.numeric=TRUE)
head(gt)
nrow(gt)
NA_prop = rowSums(is.na(gt))/ncol(gt)
plot(sort(NA_prop))

NA_prop_df = data.frame(
    index = 1:length(NA_prop),
    NA_prop = NA_prop
)
NA_prop_df$col = ifelse(NA_prop_df$NA_prop == 0, "red", "black")

no_NA_names = names(NA_prop[NA_prop == 0])
length(no_NA_names)
#1748
gt.no_na = gt[row.names(gt) %in% no_NA_names,]
nrow(gt.no_na)
gt = gt.no_na
gt[is.na(gt)] = -9
head(gt)
gt_str = t(gt)
nrow(gt_str)
#5 individuals (for structure input options)
n_snps = ncol(gt_str)
n_snps
#1748 snps  (for structure input options)
colnames(gt_str) = paste("SNP_", seq(1, n_snps), sep = "")
head(gt_str[,1:5])

#adding the loc info
sample_ids = row.names(gt_str)
metadata = read.csv("data/sample_metadata/Nc_canton_loc_date.lat_lon.csv")
loc_data = left_join(data.frame(Sequence_label = sample_ids), metadata %>% dplyr::select(Sequence_label, Canton) )
loc_int = data.frame(Canton = unique(loc_data$Canton), int = 1:length(unique(loc_data$Canton)))
loc_data.int = left_join(loc_data, loc_int, by = "Canton")
gt_str.loc = data.frame(LocData = loc_data.int$int, gt_str)
gt_str.loc[,1:5] %>% head

write.table(gt_str.loc, "data/Nc/final_tables/FINAL_snp.thinned_no_NA.LocData.structure", quote = F, sep = "\t", row.names = T, col.names = NA)
#subsequently delete the LocData colname manually or with sed
