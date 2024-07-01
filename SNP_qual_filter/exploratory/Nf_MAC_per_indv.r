library(vcfR)
library(tidyr)
#library(dplyr)
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
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)
gq <- extract.gt(vcf, element='GQ', as.numeric=TRUE)

rm(list = c("dp", "gq"))
gc()

dp.long = dp %>% data.frame %>% pivot_longer(names_to = "sample", values_to = "DP", cols = everything())
gq.long = gq %>% data.frame %>% pivot_longer(names_to = "sample", values_to = "GQ", cols = everything())

nrow(dp.long) == nrow(gq.long)
head(dp.long)
head(gq.long)

dp.long = dp.long[!is.na(gq.long$GQ),]
gq.long = gq.long[!is.na(gq.long$GQ),]
head(gq.long)
nrow(dp.long) == nrow(gq.long)

dp_gq.long = cbind(dp.long, gq.long["GQ"])
rm(list = c("dp.long", "gq.long"))
gc()


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


dp_gq.long$lib = vector(mode = "character", length = nrow(dp_gq.long))
dp_gq.long[dp_gq.long$sample %in% first_set.ids, "lib"] = "lib_1"
dp_gq.long[dp_gq.long$sample %in% second_set.ids, "lib"] = "lib_2"
dp_gq.long[dp_gq.long$sample %in% third_set.ids, "lib"] = "lib_3"
dp_gq.long[dp_gq.long$sample %in% fourth_set.ids, "lib"] = "lib_4"

dp_gq.long.subsample = sample(x = 1:nrow(dp_gq.long), size = 10^5, replace = F)

dp_gq.long = dp_gq.long[dp_gq.long.subsample,]


mod1 = lm(GQ ~ log10(DP+1) + lib, data = dp_gq.long)
qqnorm(residuals(mod1))
plot(residuals(mod1))
summary(mod1)
#DP relationship is sig, lib 3 int is sig (with no log10 of DP)
#everything is sig with log10 of DP
anova(mod1)
#P = 0.2e-9

mod2 = lm(GQ ~ log10(DP+1) * lib, data = dp_gq.long)
qqnorm(residuals(mod2))
plot(residuals(mod2))
summary(mod2)
#everything is sig...
anova(mod2)
#all sig

p1 = ggplot(
    dp_gq.long, 
    aes(x = DP, y = GQ)
) +
    facet_wrap(~lib, ncol = 1) +
    geom_point(alpha = 0.25)  +
    geom_hline(yintercept = 39.5, color = "blue") +
    geom_vline(xintercept = 2.5, color = "red") +
    scale_x_log10() +
    theme_classic() +
    labs(x = "DP", y = "GQ") +
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
