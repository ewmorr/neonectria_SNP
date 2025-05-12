library(dplyr)
library(ggplot2)
library(scales)
library(gridExtra)
source("library/ggplot_theme.txt")

#PED sample IDs
sample_ids.noNC = read.table("data/Nf/final_tables/rm_dups/FINAL_snp.gwas_analyses.cluster2_no_NC.sampleIDs", header = F)
sample_ids.wNC = read.table("data/Nf/final_tables/rm_dups/FINAL_snp.gwas_analyses.cluster2.sampleIDs", header = F)

#row_ids = which(sample_ids$V1 %in% sample_metadata.site_info$Sequence_label)

#The genotype data can simply be read in as a matrix (according the docs)
#OR can try loading LEA and using readLfmm()
Y.noNC = as.matrix(read.table("data/Nf/final_tables/rm_dups/FINAL_snp.gwas_analyses.cluster2_no_NC.lfmm", header = F))
Y.wNC = as.matrix(read.table("data/Nf/final_tables/rm_dups/FINAL_snp.gwas_analyses.cluster2.lfmm", header = F))

SNP_pos.noNC = read.table("data/Nf/final_tables/rm_dups/FINAL_snp.gwas_analyses.cluster2_no_NC.recode.map")
SNP_pos.wNC = read.table("data/Nf/final_tables/rm_dups/FINAL_snp.gwas_analyses.cluster2.recode.map")

SNP_pos.noNC = SNP_pos.noNC[c(1,4)]
SNP_pos.wNC = SNP_pos.wNC[c(1,4)]
colnames(SNP_pos.noNC) = c("scaffold", "position")
colnames(SNP_pos.wNC) = c("scaffold", "position")

scf_lens = read.table("data/Nf/final_tables/rm_dups/scaffold_lengths.csv", sep = ",", header = F)
colnames(scf_lens) = c("scaffold", "length")

SNP_pos.noNC = left_join(SNP_pos.noNC, scf_lens)
SNP_pos.wNC = left_join(SNP_pos.wNC, scf_lens)


###################
###################
###################
# missingness

########
# no NC
dim(Y.noNC)
#424811

# how many cols have missing vals
sum( colSums(Y.noNC == 9) > 0 )
sum( colSums(Y.noNC == 9) > 0 )/ncol(Y.noNC)
# 84%

range( colSums(Y.noNC == 9)/nrow(Y.noNC) )
# 0 - 31.5 percent
plot(sort( colSums(Y.noNC == 9)/nrow(Y.noNC) ))

# how many SNPs with less than 10% missing
sum( colSums(Y.noNC == 9)/nrow(Y.noNC) < 0.1 )
# 305870
# <5%
sum( colSums(Y.noNC == 9)/nrow(Y.noNC) < 0.05)
# 207693

nrow(Y.noNC)*0.05
#at 5% missing at most impute 4 samples
nrow(Y.noNC)*0.1

SNP_pos.noNC$percent_NA = colSums(Y.noNC == 9)/nrow(Y.noNC)
SNP_pos.noNC$num_NA = colSums(Y.noNC == 9)


rowSums(Y.noNC == 9)/ncol(Y.noNC)

######################################
# plotting missing ness by position

# order by size
scf_order = (scf_lens %>% filter(length > 10^5) %>% pull(scaffold))[
    order(scf_lens %>% filter(length > 10^5) %>% pull(length), decreasing = T )
]


p1 = ggplot(
    SNP_pos.noNC %>% filter(percent_NA < 0.03 & length > 10^5),
    aes(
        x = position/10^3,
        y = percent_NA*100
    )
) +
    xlim(0,NA) +
    geom_point(alpha = 0.05, size = 0.5) +
    scale_y_continuous(breaks = 0:5) +
    facet_wrap(~factor(scaffold, levels = scf_order), scales = "free_x", ncol = 4) +
    labs(x = "Position (Kb)", y = "% missing") +
    theme_bw()
p1

p2 = ggplot(
    SNP_pos.noNC %>% filter(length > 10^5),
    aes(
        x = position/10^3,
        y = percent_NA*100
    )
) +
    xlim(0,NA) +
    geom_point(alpha = 0.05, size = 0.5) +
    #scale_y_continuous(breaks = 0:5) +
    facet_wrap(~factor(scaffold, levels = scf_order), scales = "free_x", ncol = 4) +
    labs(x = "Position (Kb)", y = "% missing") +
    theme_bw()
p2

png("figures/GxE/RDA/missingness_per_scaffold.noNC.png", width = 2160, height = 1080, res = 150)
grid.arrange(p2, p1, ncol = 2)
dev.off()

nrow(SNP_pos.noNC %>% filter(length > 10^5 & percent_NA < 0.05))

p3 = ggplot(
    SNP_pos.noNC %>% filter(length > 10^5),
    aes(y = sort(percent_NA*100), x =  seq(1:length(percent_NA)))
) +
    geom_point() +
    labs(y = "% missing", x = "SNP rank") +
    theme_bw()
p3
png("figures/GxE/RDA/missingness_ranked.noNC.png", width = 1500, height = 500, res = 150)
p3
dev.off()

noNC_keep_SNPs = which(SNP_pos.noNC$percent_NA < 0.05 & SNP_pos.wNC$length > 10^5)

#########################################
# here is the filtered DF
Y.noNC.lt5perc = Y.noNC[,noNC_keep_SNPs]
dim(Y.noNC.lt5perc)
sort(rowSums(Y.noNC.lt5perc == 9)/ncol(Y.noNC.lt5perc))
sample_ids.noNC[which(rowSums(Y.noNC.lt5perc == 9)/ncol(Y.noNC.lt5perc)>0.05),]
#NH.bart has two
#CT has 2
# otherwise all are single site
# this distribution to much with the stats too much

SNP_pos.noNC.lt5perc = SNP_pos.noNC[noNC_keep_SNPs,]
nrow(SNP_pos.noNC.lt5perc)

#########################################

########
# with NC
dim(Y.wNC)
#424811

# how many cols have missing vals
sum( colSums(Y.wNC == 9) > 0 )
sum( colSums(Y.wNC == 9) > 0 )/ncol(Y.wNC)
# 85%

range( colSums(Y.wNC == 9)/nrow(Y.wNC) )
# 0 - 29.6 percent
plot(sort( colSums(Y.wNC == 9)/nrow(Y.wNC) ))

# how many SNPs with less than 10% missing
sum( colSums(Y.wNC == 9)/nrow(Y.wNC) < 0.1 )
# 290569
# <5%
sum( colSums(Y.wNC == 9)/nrow(Y.wNC) < 0.05)
# 195223

nrow(Y.wNC)*0.05
#at 5% missing at most impute 4 samples
nrow(Y.wNC)*0.1

SNP_pos.wNC$percent_NA = colSums(Y.wNC == 9)/nrow(Y.wNC)
SNP_pos.wNC$num_NA = colSums(Y.wNC == 9)


sort(rowSums(Y.wNC == 9)/ncol(Y.wNC))

######################################
# plotting missing ness by position

# order by size
scf_order = (scf_lens %>% filter(length > 10^5) %>% pull(scaffold))[
    order(scf_lens %>% filter(length > 10^5) %>% pull(length), decreasing = T )
]


p1 = ggplot(
    SNP_pos.wNC %>% filter(percent_NA < 0.05 & length > 10^5),
    aes(
        x = position/10^3,
        y = percent_NA*100
    )
) +
    xlim(0,NA) +
    geom_point(alpha = 0.05, size = 0.5) +
    scale_y_continuous(breaks = 0:5) +
    facet_wrap(~factor(scaffold, levels = scf_order), scales = "free_x", ncol = 4) +
    labs(x = "Position (Kb)", y = "% missing") +
    theme_bw()
p1

p2 = ggplot(
    SNP_pos.wNC %>% filter(length > 10^5),
    aes(
        x = position/10^3,
        y = percent_NA*100
    )
) +
    xlim(0,NA) +
    geom_point(alpha = 0.05, size = 0.5) +
    #scale_y_continuous(breaks = 0:5) +
    facet_wrap(~factor(scaffold, levels = scf_order), scales = "free_x", ncol = 4) +
    labs(x = "Position (Kb)", y = "% missing") +
    theme_bw()
p2

png("figures/GxE/RDA/missingness_per_scaffold.wNC.png", width = 2160, height = 1080, res = 150)
grid.arrange(p2, p1, ncol = 2)
dev.off()

nrow(SNP_pos.wNC %>% filter(length > 10^5 & percent_NA < 0.05))

p3 = ggplot(
    SNP_pos.wNC %>% filter(length > 10^5),
    aes(y = sort(percent_NA*100), x =  seq(1:length(percent_NA)))
) +
    geom_point() +
    labs(y = "% missing", x = "SNP rank") +
    theme_bw()
p3
png("figures/GxE/RDA/missingness_ranked.wNC.png", width = 1500, height = 500, res = 150)
p3
dev.off()

wNC_keep_SNPs = which(SNP_pos.wNC$percent_NA < 0.05 & SNP_pos.wNC$length > 10^5)

#########################################
# here is the filtered DF
Y.wNC.lt5perc = Y.wNC[,wNC_keep_SNPs]
dim(Y.wNC.lt5perc)
sort(rowSums(Y.wNC.lt5perc == 9)/ncol(Y.wNC.lt5perc))
sample_ids.wNC[which(rowSums(Y.wNC.lt5perc == 9)/ncol(Y.wNC.lt5perc)>0.05),]
#NH.bart has two
#CT has 2
# otherwise all are single site
# this distribution to much with the stats too much

SNP_pos.wNC.lt5perc = SNP_pos.wNC[wNC_keep_SNPs,]
nrow(SNP_pos.wNC.lt5perc)

#########################################

#########################################
# IMPUTATION BASED ON MAX VALUE PER SNP
#overall missingness for each
sum(Y.noNC.lt5perc == 9)/(nrow(Y.noNC.lt5perc) * ncol(Y.noNC.lt5perc))
# overall missingness 0.01695916
# max 5% (4 samples) per SNP
sum(Y.wNC.lt5perc == 9)/(nrow(Y.wNC.lt5perc) * ncol(Y.wNC.lt5perc))
# overall missingness 0.01601603
# max 5% (4 samples) per SNP
dim(Y.noNC.lt5perc)

# convert val == 9 to NA
Y.noNC.lt5perc.impute = Y.noNC.lt5perc
Y.noNC.lt5perc.impute[Y.noNC.lt5perc.impute == 9] = NA
sum(is.na(Y.noNC.lt5perc.impute))

Y.noNC.lt5perc.impute <- apply(
    Y.noNC.lt5perc.impute, 
    2, 
    function(x) replace(
        x, 
        is.na(x), 
        as.numeric(names(which.max(table(x))))) # tables the counts of 0/1 and takes the label of the most frequent then conversts back to numeric
)
sum(is.na(Y.noNC.lt5perc.impute))

Y.wNC.lt5perc.impute = Y.wNC.lt5perc
Y.wNC.lt5perc.impute[Y.wNC.lt5perc.impute == 9] = NA
sum(is.na(Y.wNC.lt5perc.impute))

Y.wNC.lt5perc.impute <- apply(
    Y.wNC.lt5perc.impute, 
    2, 
    function(x) replace(
        x, 
        is.na(x), 
        as.numeric(names(which.max(table(x))))) # tables the counts of 0/1 and takes the label of the most frequent then converts back to numeric
)
sum(is.na(Y.wNC.lt5perc.impute))

#########################################
dim(Y.noNC.lt5perc.impute)
dim(Y.wNC.lt5perc.impute)
#########################################

#########################################
#MAC filter > 0.05

nrow(Y.noNC.lt5perc.impute) * 0.05
nrow(Y.wNC.lt5perc.impute) * 0.05
# gt 4 samples in both cases
nrow(SNP_pos.noNC.lt5perc) == ncol(Y.noNC.lt5perc.impute)
nrow(SNP_pos.wNC.lt5perc) == ncol(Y.wNC.lt5perc.impute)

#filter for MAC ge 4
minMAC = 4

###################
# no NC
ref_sum = colSums(Y.noNC.lt5perc.impute == 0)
alt_sum = colSums(Y.noNC.lt5perc.impute == 1)

nrow(Y.noNC.lt5perc.impute)
sum(ref_sum < minMAC)
# 88698
sum(alt_sum < minMAC)
ncol(Y.noNC.lt5perc.impute)
# 207335
ncol(Y.noNC.lt5perc.impute) - sum(ref_sum < minMAC)
# 118637

which(ref_sum < minMAC) 
which(alt_sum < minMAC)# there are none
rm_cols = which(ref_sum < minMAC)

length(rm_cols)
#88698
length(rm_cols)/ncol(Y.noNC.lt5perc.impute)
# 0.4278004
# 
# ##########################################
#MAC filtered
Y.noNC.filteredMAC = Y.noNC.lt5perc.impute[,-rm_cols]
SNP_pos.noNC.filteredMAC = SNP_pos.noNC.lt5perc[-rm_cols,]
ncol(Y.noNC.filteredMAC)
nrow(SNP_pos.noNC.filteredMAC)


p1 = ggplot(
    SNP_pos.noNC.filteredMAC %>% 
    group_by(scaffold) %>%
    arrange(position) %>%
    mutate(diff = position - lag(position, default = first(0)) )
,
    aes(
        x = position/10^3,
        y = diff
    )
) +
    xlim(0,NA) +
    geom_point(alpha = 0.15, size = 0.5) +
    scale_y_log10(breaks = c(1,10,10^2,10^3,10^4,10^5)) +
    facet_wrap(~factor(scaffold, levels = scf_order), scales = "free", ncol = 2) +
    labs(x = "Position (Kb)", y = "distance n:n-1 (bb)") +
    theme_bw()
p1

png("figures/GxE/RDA/SNP_distance.noNC.NAlt5perc_and_MACgt5perc.png", width = 2160, height = 1400, res = 150)
p1
dev.off()


###################
# w NC
ref_sum = colSums(Y.wNC.lt5perc.impute == 0)
alt_sum = colSums(Y.wNC.lt5perc.impute == 1)

nrow(Y.wNC.lt5perc.impute)
sum(ref_sum < minMAC)
# 81222
sum(alt_sum < minMAC)
ncol(Y.wNC.lt5perc.impute)
# 194873
ncol(Y.wNC.lt5perc.impute) - sum(ref_sum < minMAC)
# 113651

which(ref_sum < minMAC) 
which(alt_sum < minMAC)# there are none
rm_cols = which(ref_sum < minMAC)

length(rm_cols)
#81222
length(rm_cols)/ncol(Y.wNC.lt5perc.impute)
# 0.4167945
# 
# ##########################################
#MAC filtered
Y.wNC.filteredMAC = Y.wNC.lt5perc.impute[,-rm_cols]
SNP_pos.wNC.filteredMAC = SNP_pos.wNC.lt5perc[-rm_cols,]
ncol(Y.wNC.filteredMAC)
nrow(SNP_pos.wNC.filteredMAC)


p1 = ggplot(
    SNP_pos.wNC.filteredMAC %>% 
    group_by(scaffold) %>%
    arrange(position) %>%
    mutate(diff = position - lag(position, default = first(0)) )
,
    aes(
        x = position/10^3,
        y = diff
    )
) +
    xlim(0,NA) +
    geom_point(alpha = 0.15, size = 0.5) +
    scale_y_log10(breaks = c(1,10,10^2,10^3,10^4,10^5)) +
    facet_wrap(~factor(scaffold, levels = scf_order), scales = "free", ncol = 2) +
    labs(x = "Position (Kb)", y = "distance n:n-1 (bp)") +
    theme_bw()
p1

png("figures/GxE/RDA/SNP_distance.wNC.NAlt5perc_and_MACgt5perc.png", width = 2160, height = 1400, res = 150)
p1
dev.off()

###############################
# write filtered files for RDA
###############################

write.table(Y.noNC.filteredMAC, "data/Nf/GxE/RDA/Y.noNC.filteredMAC.tsv", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(Y.wNC.filteredMAC, "data/Nf/GxE/RDA/Y.wNC.filteredMAC.tsv", row.names = F, col.names = F, quote = F, sep = "\t")
write.csv(SNP_pos.noNC.filteredMAC, "data/Nf/GxE/RDA/SNP_pos.noNC.filteredMAC.csv", row.names = F, quote = F)
write.csv(SNP_pos.wNC.filteredMAC, "data/Nf/GxE/RDA/SNP_pos.wNC.filteredMAC.csv", row.names = F, quote = F)
write.table(sample_ids.noNC, "data/Nf/GxE/RDA/sample_ids.noNC", row.names = F, col.names = F, quote = F)
write.table(sample_ids.wNC, "data/Nf/GxE/RDA/sample_ids.wNC", row.names = F, col.names = F, quote = F)

