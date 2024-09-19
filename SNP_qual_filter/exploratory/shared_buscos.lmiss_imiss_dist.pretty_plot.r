library(ggplot2)
library(gridExtra)

lmiss = read.table("data/shared_buscos/INFOfilters.removed.DPGQ_filter.lmiss", header = T)
imiss = read.table("data/shared_buscos/INFOfilters.removed.DPGQ_filter.imiss", header = T)

head(lmiss)
head(imiss)

p1 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count (SNP)") 
p1

p1 = p1 + 
    geom_vline(xintercept = 0.25*100, color = "red") +
    expand_limits(x = c(0,100), y = c(0,5.5e+05))  +
    annotate(geom = "text", x = (0.25*100)+1, hjust = 0, y = 3e+05, label = "25%") +
    labs(title = "Before filtering for missingness")
p1

p2 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)", y = "count (indv)", title = "") +
    expand_limits(x = c(0,30), y = c(0,30))  +
    scale_x_continuous(breaks = c(0,15,30,45,60)) +
    labs(title = "")
p2


###########################################
#First pass -- after remove top 5% loci -- Import new tables and rerun above

lmiss = read.table("data/shared_buscos/NA_filter.loc_gt_0.25.lmiss", header = T)
imiss = read.table("data/shared_buscos/NA_filter.loc_gt_0.25.imiss", header = T)

p3 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count (SNP)") +
    expand_limits(x = c(0,100), y = c(0,5.5e+05))  +
    labs(title = "After locus 25% percentile filter")
p3

p4 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)", y = "count (indv)", title = "") +
    expand_limits(x = c(0,30), y = c(0,30))  +
    scale_x_continuous(breaks = c(0,15,30,45,60)) 

p4




pdf("figures/quality_filtering/shared_buscos/shared_buscos.lmiss_imiss_distribution.locus_hard_filter.pdf", height = 4, width = 8)
grid.arrange(p1,p2,p3,p4, ncol = 2, widths = c(0.5,0.5))
dev.off()

