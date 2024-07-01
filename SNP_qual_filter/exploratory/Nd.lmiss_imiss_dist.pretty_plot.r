library(ggplot2)
library(gridExtra)

lmiss = read.table("data/Nd/all_libs.DPGQ_filter.lmiss", header = T)
imiss = read.table("data/Nd/all_libs.DPGQ_filter.imiss", header = T)

head(lmiss)
head(imiss)

p1 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 30) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count (SNP)") +
    expand_limits(x = c(0,100), y = c(0,6e+05))  +
    labs(title = "Before filtering for missingness")
p1


p2 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)", y = "count (indv)", title = "") +
    expand_limits(x = c(0,60), y = c(0,30))  +
    scale_x_continuous(breaks = c(0,15,30,45,60)) +
    geom_vline(xintercept = 0.15*100, color = "red") +
    annotate(geom = "text", x = (0.15*100)+1, hjust = 0, y = 10, label = ">15%") +
    labs(title = "")
p2


###########################################
#First pass -- there is one individual that is above 60% missingness after 
# applying the DP and GQ filters, whereas all others are below 15%
# let's remove this (so everything now below 15% like Nf) and then apply a 
# single hard filter to the SNPs to 25%

lmiss = read.table("data/Nd/NA_filter.ind_gt_0.15.lmiss", header = T)
imiss = read.table("data/Nd/NA_filter.ind_gt_0.15.imiss", header = T)

p3 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins =30) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count (SNP)") +
    #geom_vline(xintercept = quantile(lmiss$F_MISS, c(0.95))*100, color = "red") +
    #annotate(geom = "text", x = (quantile(lmiss$F_MISS, c(0.95))*100)+1, hjust = 0, y = 3e+05, label = "95th %") +
    geom_vline(xintercept = 0.25*100, color = "red") +
    annotate(geom = "text", x = (0.25*100)+1, hjust = 0, y = 3e+05, label = ">25%") +
    expand_limits(x = c(0,100), y = c(0,6e+05))  +
    labs(title = "After 1st indv 15% filter")
p3

p4 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)", y = "count (indv)", title = "") +
    expand_limits(x = c(0,60), y = c(0,30))  +
    scale_x_continuous(breaks = c(0,15,30,45,60)) 

p4


#############################################
#Second pass -- after removed sites >25% missing


lmiss = read.table("data/Nd/NA_filter.loc_gt_0.25.ind_gt_0.15.lmiss", header = T)
imiss = read.table("data/Nd/NA_filter.loc_gt_0.25.ind_gt_0.15.imiss", header = T)

p5 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins =30) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count (SNP)") +
    expand_limits(x = c(0,100), y = c(0,6e+05))  +
    #geom_vline(xintercept = quantile(lmiss$F_MISS, c(0.95))*100, color = "red") +
    #annotate(geom = "text", x = (quantile(lmiss$F_MISS, c(0.95))*100)+1, hjust = 0, y = 3e+05, label = "95th %") +
    labs(title = "Final loci <25% + indv <15%")
p5

p6 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)", y = "count (indv)", title = "") +
    scale_x_continuous(breaks = c(0,15,30,45,60)) +
    expand_limits(x = c(0,60), y = c(0,30)) 
p6


pdf("figures/quality_filtering/Nd/Nd.lmiss_imiss_distribution.hard_filter.pdf", height = 5.5, width = 8)
grid.arrange(p1,p2,p3,p4,p5,p6, ncol = 2, widths = c(0.5,0.5))
dev.off()

