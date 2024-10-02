library(vcfR)

# Nf
vcf <- read.vcfR("data/Nf/final_tables/rm_dups/FINAL_snp.mac_ge2.biallele.LD.structure_analyses.vcf.gz", verbose = FALSE)
gt <- extract.gt(vcf, element='GT', as.numeric=TRUE)
head(gt)
gt[is.na(gt)] = -9
head(gt)
gt_str = t(gt)
nrow(gt_str)
#115 individuals (for structure input options)
n_snps = ncol(gt_str)
n_snps
#188237 snps  (for structure input options)
colnames(gt_str) = paste("SNP_", seq(1, n_snps), sep = "")
head(gt_str[,1:5])
write.table(gt_str, "data/Nf/final_tables/rm_dups/FINAL_snp.structure", quote = F, sep = "\t", row.names = T, col.names = NA)

#thinned Nf
vcf <- read.vcfR("data/Nf/final_tables/rm_dups/FINAL_snp.mac_ge2.biallele.LD.structure_analyses.vcf.gz", verbose = FALSE)
gt <- extract.gt(vcf, element='GT', as.numeric=TRUE)
head(gt)
nrow(gt)
NA_prop = rowSums(is.na(gt))/ncol(gt)
plot(sort(NA_prop))
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
write.table(gt_str, "data/Nf/final_tables/rm_dups/FINAL_snp.thinned_no_NA.structure", quote = F, sep = "\t", row.names = T, col.names = NA)


# Nd
vcf <- read.vcfR("data/Nd/final_tables/rm_dups/FINAL_snp.mac_ge2.biallele.LD.structure_analyses.vcf.gz", verbose = FALSE)
gt <- extract.gt(vcf, element='GT', as.numeric=TRUE)
head(gt)
gt[is.na(gt)] = -9
head(gt)
gt_str = t(gt)
nrow(gt_str)
#30 individuals (for structure input options)
n_snps = ncol(gt_str)
n_snps
#81334 snps  (for structure input options)
colnames(gt_str) = paste("SNP_", seq(1, n_snps), sep = "")
head(gt_str[,1:5])
write.table(gt_str, "data/Nd/final_tables/rm_dups/FINAL_snp.structure", quote = F, sep = "\t", row.names = T, col.names = NA)



# Nc
vcf <- read.vcfR("data/Nc/final_tables/FINAL_snp.mac_ge2.biallele.LD.structure_analyses.vcf.gz", verbose = FALSE)
gt <- extract.gt(vcf, element='GT', as.numeric=TRUE)
head(gt)
gt[is.na(gt)] = -9
head(gt)
gt_str = t(gt)
nrow(gt_str)
#5 individuals (for structure input options)
n_snps = ncol(gt_str)
n_snps
#2210 snps  (for structure input options)
colnames(gt_str) = paste("SNP_", seq(1, n_snps), sep = "")
head(gt_str[,1:5])
write.table(gt_str, "data/Nc/final_tables/FINAL_snp.structure", quote = F, sep = "\t", row.names = T, col.names = NA)
