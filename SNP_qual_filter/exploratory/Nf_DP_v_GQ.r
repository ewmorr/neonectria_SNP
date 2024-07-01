library(vcfR)
library(tidyr)
library(ggplot2)
source("ggplot_theme.txt")

#filtered VCF
vcf <- read.vcfR("data/Nf/INFOfilters.removed.recode.vcf", verbose = FALSE)

####################
#Site level metrics#
#i.e., locus x indv #
####################

#calculate MAC per individual and test for lib effects
#dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)
gt <- extract.gt(vcf, element='GT', as.numeric=TRUE)
head(gt)

gt[rowSums(gt, na.rm = T) == 0,] %>% nrow()
#339385
#there are some invariant sites after removing the one indv
# so will need to add this to any calcs from vcftools which picks up 342658
# this does not agree with our calc above... there is likely something funny 
# going on with the way vcftools is handling NA values, 
# e.g., if there is only one SNP called (and rest NA) per row is it a MAC=1?

gt[rowSums(gt, na.rm = T) == 1 & rowSums(is.na(gt)) != 123,] %>% nrow()
#290639 
gt[rowSums(gt, na.rm = T) == 123 & rowSums(is.na(gt)) == 0,] %>% nrow()
#1551
gt[rowSums(gt, na.rm = T) == 123 & rowSums(is.na(gt)) == 1,] %>% nrow()
#86 
#Oh... there are some oddball cases where all but one are ALT not REF. 
# This is probably bc of removing the one sample.
#292441 total so still doesn't add up to vcftools
gt[rowSums(gt, na.rm = T) == 0 & rowSums(is.na(gt)) == 123,] %>% nrow()
#0
# there are no calls with one REF and rest NA
gt[rowSums(is.na(gt)) == 123,] %>% nrow()
#127
#there are 127 rows total with all but one NA
gt[rowSums(is.na(gt)) == 123 & rowSums(gt, na.rm = T) == 1,] %>% rowSums(.,na.rm = T) %>% sum()
#127
#all of the 123 NA samples are ALTs (interesting...)
gt[rowSums(is.na(gt)) == 124,] %>% nrow()
#2 with all NAs
gt[rowSums(is.na(gt)) == 0,] %>% nrow()
#684333 with no NA
gt[rowSums(gt, na.rm = T) == 0,] %>% nrow()
#339385 are all REF
gt[rowSums(gt, na.rm = T) == 124,] %>% nrow()
#518 are all ALT
gt[rowSums(gt, na.rm = T) + rowSums(is.na(gt)) ==  124,] %>% nrow()
#973 are all NA and ALT

973 + 518 + 339385 + 2 + 127 + 1551 + 86 + 290639
# 633281
# above is getting closer to the number of SNPs removed with vcftools --mac 2
# filter, which removes 675684. 

gt[rowSums(gt, na.rm = T) + rowSums(is.na(gt)) >  124,] %>% nrow()
#3932 are > 2 alleles... (ALT beyond biallele are coded sequentially...)
633281 + 3932
#637213

#well...in any case we are not getting the same numbers as vcftools; 
# there is likely a category we are missing above or a funky calculation in vcftools
# but doubtful we are biasing the below MAC == 1 comparison

nrow(gt)
#1322180

#filter to only SNPs with mac == 1
#first get the rows where REF is the minor allele and invert the values 
#to make the following calcs easier
gt.altmac1 = gt[rowSums(gt, na.rm = T) == 123 & rowSums(is.na(gt)) == 0,]
gt.altmac1[gt.altmac1 == 1] = -9
gt.altmac1[gt.altmac1 == 0] = 1
gt.altmac1[gt.altmac1 == -9] = 0
nrow(gt.altmac1)
rowSums(gt.altmac1) %>% sum()
gt.altmac1[rowSums(gt.altmac1, na.rm = T) != 1 ,] %>% nrow()
# 27 rows have multialleles that are messing up the calcs
gt.altmac1[rowSums(gt.altmac1 <= 1) != 124,] %>% nrow()
#remove them
gt.altmac1 = gt.altmac1[rowSums(gt.altmac1, na.rm = T) == 1,]


gt.mac1 = rbind(gt[rowSums(gt, na.rm = T) == 1 & rowSums(is.na(gt)) != 123,],  gt.altmac1)
ncol(gt.mac1)
nrow(gt.mac1)

head(gt.mac1)
#sum by individual to get counts per indv
indv.mac1 = colSums(gt.mac1, na.rm = T)
head(indv.mac1)

#sample_ids
first_set.ids = paste0("NG", 1:99)
second_set.ids = paste0("NG", 101:163)
third_set.ids = c(
"NG103",
"NG111",
"NG112",
"NG114",
"NG116",
"NG121",
"NG155",
"NG160",
"NG117",
"NG163"
)
second_set.ids = second_set.ids[!second_set.ids %in% third_set.ids]
fourth_set.ids = paste0("NG", 170:196)



ind_stats = data.frame(
    ind_names = names(indv.mac1),
    ind_mac1 = indv.mac1
    )
ind_stats$lib = vector(mode = "character", length = nrow(ind_stats))
ind_stats[row.names(ind_stats) %in% first_set.ids, "lib"] = "lib_1"
ind_stats[row.names(ind_stats) %in% second_set.ids, "lib"] = "lib_2"
ind_stats[row.names(ind_stats) %in% third_set.ids, "lib"] = "lib_3"
ind_stats[row.names(ind_stats) %in% fourth_set.ids, "lib"] = "lib_4"

mod1 = aov(ind_mac1 ~ lib, data = ind_stats)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
#P = 0.167
TukeyHSD(mod1)
#no difs


p1 = ggplot(
    ind_stats, 
    aes(x = reorder(ind_names, ind_mac1), y = ind_mac1)
) +
    geom_col() +
    theme_classic() +
    labs(x = "Individual", y = "MAC == 1 SNP count") +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, size = 6)
    )

p1

p2 = ggplot(
    ind_stats, 
    aes(x = reorder(ind_names, ind_mac1), y = ind_mac1)
) +
    geom_col() +
    theme_classic() +
    labs(x = "Individual", y = "MAC == 1 SNP count") +
    facet_wrap(~lib, ncol = 1) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, size = 6)
    )

p2

p3 = ggplot(
    ind_stats, 
    aes(x = ind_mac1, fill = lib)
) +
    geom_histogram(bins = 20) +
    theme_classic()  +
    scale_fill_manual(values = cbPalette) +
    labs(x = "Singleton alleles per individual", fill = "Library", y = "Number individuals")

p3

pdf("figures/quality_filtering/Nf/MAC_singletons_per_indv_by_library.pdf", width = 12, height = 6)
p1
p2
p3
dev.off()
