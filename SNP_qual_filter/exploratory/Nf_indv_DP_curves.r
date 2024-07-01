library(vcfR)
library(tidyr)
library(ggplot2)
source("ggplot_theme.txt")

#filtered VCF
vcf <- read.vcfR("data/Nf/out.vcf", verbose = FALSE)

####################
#Site level metrics#
#i.e., locus x indv #
####################

#after filtering for variant-level quality, filter at the site level based on 
# DP and missingness
#plots of initial mean DP and percent missing 

dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)
gt <- extract.gt(vcf, element='GT', as.numeric=TRUE)
dp[is.na(gt)] = NA
rm(gt)
rm(vcf)
gc() # "garbage collection" to dump the mem
#########################
# there is a bimodal distribution in locus mean DP
# but only in libs one and two
# extract those samples and dump the rest to save mem

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

dp.one = dp[,colnames(dp) %in% first_set.ids]
dp.two = dp[,colnames(dp) %in% second_set.ids]

rm(dp)
gc()



#plot all depth curves by individual

dp.one.long = pivot_longer(
    data.frame(genotype = rownames(dp.one), dp.one), 
    cols = -genotype, names_to = "indv", values_to = "DP"
)

dp.one.long[dp.one.long == 0] = NA
dp.one.long = na.omit(dp.one.long)

dp.two.long = pivot_longer(
    data.frame(genotype = rownames(dp.two), dp.two), 
    cols = -genotype, names_to = "indv", values_to = "DP"
)

dp.two.long[dp.two.long == 0] = NA
dp.two.long = na.omit(dp.two.long)

rm(list = c("dp.one", "dp.two"))
gc()

p1 = ggplot(dp.one.long, aes(DP)) +
    geom_bar() +
    theme_classic() +
    labs(x = "Locus DP\n(values of above 100 excluded)") +
    facet_wrap(~indv, scales = "free_y") +
    scale_x_continuous(limits = c(0,100), breaks = c(0,25,50,50,75,100)) +
    geom_vline(xintercept = 6) + #approcimate peaks in the library one locus means graphs
    geom_vline(xintercept = 12)

pdf("figures/quality_filtering/Nf/individual_DP_curves.lib_one.pdf", width = 16, height = 12)
p1
dev.off()

rm(list = c("dp.one.long", "p1"))
gc()

p1 = ggplot(dp.two.long, aes(DP)) +
    geom_bar() +
    theme_classic() +
    labs(x = "Locus DP\n(values of above 100 excluded)") +
    facet_wrap(~indv, scales = "free_y") +
    scale_x_continuous(limits = c(0,100), breaks = c(0,25,50,50,75,100)) +
    geom_vline(xintercept = 22) + #approcimate peaks in the library two locus means graphs
    geom_vline(xintercept = 30)

pdf("figures/quality_filtering/Nf/individual_DP_curves.lib_two.pdf", width = 16, height = 12)
p1
dev.off()

rm(list = c("dp.two.long", "p1"))
gc()

# The graphs confirm that there are two or more different sample medians clustered around the lib-wide peaks, hence the bimodal distribution in the libs overall. There are two samples (NG113 and NG154) in lib two that still appear to have bimodal dist which might suggest these are dikaryons/diploid cultures. Let's keep an eye to see if these get filtered out based on NA filtering as they should have a high proportion of NAs due to excluding het calls.