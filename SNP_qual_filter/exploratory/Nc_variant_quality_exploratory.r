library(vcfR)
library(tidyr)
library(ggplot2)
source("library/ggplot_theme.txt")

#filtered VCF
vcf <- read.vcfR("data/Nc/out.vcf", verbose = FALSE)

##########################
#Variant level metrics   #
#i.e., across all samples#
##########################
#pull VCF INFO 
#alternatively split this off using cut and read in separately
#command takes a while and would save mem if reading in only info field
vcf.info = INFO2df(vcf)
vcf.qual = getQUAL(vcf)
head(vcf.qual)
vcf.qual = data.frame(QUAL = vcf.qual)

#QD
ggplot(vcf.info, aes(x = QD)) +
    geom_density() +
    theme_bw()

p1 = ggplot(vcf.info, aes(x = QD)) +
    geom_histogram() +
    theme_bw() +
    geom_vline(xintercept = 2)

#MQ
ggplot(vcf.info, aes(x = MQ)) +
    geom_density() +
    theme_bw()

p2 = ggplot(vcf.info, aes(x = MQ)) +
    geom_histogram() +
    geom_vline(xintercept = 40) +
    theme_bw()

#FS
ggplot(vcf.info, aes(x = FS)) +
    geom_density() +
    theme_bw()

p3 = ggplot(vcf.info, aes(x = FS)) +
    geom_histogram() +
    scale_x_continuous(limits = c(-2,60)) +
    theme_bw()

#DP
ggplot(vcf.info, aes(x = DP)) +
    geom_density() +
    theme_bw()

p4 = ggplot(vcf.info, aes(x = DP)) +
    geom_histogram(binwidth = 5) +
    scale_x_continuous(limits = c(0,500)) +
    theme_bw()
p4
#QUAL
p5 = ggplot(vcf.qual, aes(x = QUAL)) +
    geom_histogram() +
    #scale_x_continuous(limits = c(0,7500)) +
    theme_bw()
p5
p6 = ggplot(vcf.qual, aes(x = QUAL)) +
    geom_histogram() +
    scale_x_continuous(limits = c(0,100)) +
    theme_bw() +
    labs(x = "QUAL (>100 excluded")
p6

pdf("figures/quality_filtering/Nc/genotype_qual.SNP.pdf", width = 6, height = 4)
p1
p2
p3
p4
p5
p6
dev.off()


####################
#Site level metrics#
#i.e., locus x indv #
####################

#after filtering for variant-level quality, filter at the site level based on 
# DP and missingness
#plots of initial mean DP and percent missing 

dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)
gt <- extract.gt(vcf, element='GT', as.numeric=TRUE)
is.matrix(dp)


# percent missing on locus and ind
locus_na = (is.na(gt) %>% rowSums())/ncol(gt)
ind_na = (is.na(gt) %>% colSums())/nrow(gt)

#also take means of dp at both geno and ind to get target level for filtering
# individual sites (above we took total DP at locus, i.e., INFO field)
locus_dp_mean = rowSums(dp, na.rm = T)/ncol(dp)
ind_dp_mean = colSums(dp, na.rm = T)/nrow(dp)


ind_stats = data.frame(
    ind_names = names(ind_na),
    ind_na = ind_na,
    ind_dp_mean = ind_dp_mean
    )
rownames(ind_stats) = ind_stats$ind_names

p1 = ggplot(
    data.frame(locus_na = locus_na), 
    aes(x = locus_na*100)
) +
    geom_bar() +
    theme_classic() +
    labs(x = "Locus % missing")

p1

p2 = ggplot(
    ind_stats, 
    aes(y = ind_na*100, x = reorder(ind_names, ind_na))
) +
    geom_col() +
    theme_classic() +
    labs( y = "Individual % missing") +
    theme(
        axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank()
    )
p2

range(locus_dp_mean)
median(locus_dp_mean)
mean(locus_dp_mean)
sd(locus_dp_mean)
mean(locus_dp_mean) + (2*sd(locus_dp_mean))
mean(locus_dp_mean) - (1*sd(locus_dp_mean))
locus_dp_mean[order(locus_dp_mean, decreasing = T)] 

p3 = ggplot(
    data.frame(locus_dp_mean = locus_dp_mean), 
    aes(x = locus_dp_mean)
) +
    geom_bar() +
    theme_classic() +
    labs(x = "Locus mean DP\n(values of above 100 excluded)") +
    scale_x_continuous(limits = c(0,100), breaks = c(0,5,10,15,20,25,50,50,75,100))
p3

p4 = ggplot(
    ind_stats, 
    aes(y = ind_dp_mean, x = reorder(ind_names, ind_dp_mean))
) +
    geom_col() +
    theme_classic() +
    labs(y = "Individual mean DP") +
    theme(
        axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank()
    )
p4


pdf("figures/quality_filtering/Nc/locus.NA_DP.pdf", width = 10, height = 6)
p1
p3
dev.off()

pdf("figures/quality_filtering/Nc/individual.NA_DP.pdf", width = 16, height = 6)
p2
p4
dev.off()

#########################
# there is a bimodal distribution in locus mean DP (as there is for Nf)
# Check for individual samples

dp.long = pivot_longer(data.frame(dp), names_to = "sample", values_to = "DP", cols = everything())
head(dp.long)

p5 = ggplot(dp.long, aes(x = DP)) +
    geom_histogram(binwidth = 5) +
    scale_x_continuous(limits = c(0,100)) + #there are a couple samples with a bit over 100 but not worth worrying about
    facet_wrap(~sample, ncol = 3) +
    my_gg_theme
               
pdf("figures/quality_filtering/Nc/individual_DP.facets.pdf", width = 8, height = 6)
p5
dev.off()
#no bimodal dists when looking at indiviodual samples
