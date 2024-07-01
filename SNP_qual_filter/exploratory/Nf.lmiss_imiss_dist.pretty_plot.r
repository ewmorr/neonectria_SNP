library(ggplot2)
library(gridExtra)

lmiss = read.table("data/Nf/filtering_intermediates/all_libs.DPGQ_filter.lmiss", header = T)
imiss = read.table("data/Nf/filtering_intermediates/all_libs.DPGQ_filter.imiss", header = T)

head(lmiss)
head(imiss)

p1 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count (SNP)") +
    expand_limits(x = c(0,100), y = c(0,6e+05)) 
p1

p1 = p1 + 
    geom_vline(xintercept = quantile(lmiss$F_MISS, c(0.95))*100, color = "red") +
    annotate(geom = "text", x = (quantile(lmiss$F_MISS, c(0.95))*100)+1, hjust = 0, y = 3e+05, label = "95th %") +
    labs(title = "Before filtering for missingness")

p2 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)", y = "count (indv)", title = "") +
    expand_limits(x = c(0,60), y = c(0,30))  +
    scale_x_continuous(breaks = c(0,15,30,45,60)) +
    labs(title = "")
p2


###########################################
#First pass -- after remove top 5% loci -- Import new tables and rerun above

lmiss = read.table("data/Nf/filtering_intermediates/NA_filter.loc_gt_0.524194.lmiss", header = T)
imiss = read.table("data/Nf/filtering_intermediates/NA_filter.loc_gt_0.524194.imiss", header = T)

p3 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count (SNP)") +
    expand_limits(x = c(0,100), y = c(0,6e+05))  +
    labs(title = "After 1st locus 95th percentile filter")
p3

p4 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)", y = "count (indv)", title = "") +
    expand_limits(x = c(0,60), y = c(0,30))  +
    scale_x_continuous(breaks = c(0,15,30,45,60)) +
    geom_vline(xintercept = quantile(imiss$F_MISS, c(0.99))*100, color = "red") +
    annotate(geom = "text", x = (quantile(imiss$F_MISS, c(0.99))*100)+1, hjust = 0, y = 20, label = "99th %")

p4


#############################################
#Second pass -- after removed top indv NA (NG 139)


lmiss = read.table("data/Nf/filtering_intermediates/NA_filter.loc_gt_0.524194.ind_gt_0.30502884.lmiss", header = T)
imiss = read.table("data/Nf/filtering_intermediates/NA_filter.loc_gt_0.524194.ind_gt_0.30502884.imiss", header = T)

p5 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count (SNP)") +
    expand_limits(x = c(0,100), y = c(0,6e+05))  +
    geom_vline(xintercept = quantile(lmiss$F_MISS, c(0.95))*100, color = "red") +
    annotate(geom = "text", x = (quantile(lmiss$F_MISS, c(0.95))*100)+1, hjust = 0, y = 3e+05, label = "95th %") +
    labs(title = "After 1st locus filter + 1st indv 99th percentile filter")
p5

p6 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)", y = "count (indv)", title = "") +
    scale_x_continuous(breaks = c(0,15,30,45,60)) +
    expand_limits(x = c(0,60), y = c(0,30)) 
p6

#############################################
#Third pass -- after remove top 5% loci

lmiss = read.table("data/Nf/filtering_intermediates/NA_filter.loc_gt_0.3606560.ind_gt_0.30502884.lmiss", header = T)
imiss = read.table("data/Nf/filtering_intermediates/NA_filter.loc_gt_0.3606560.ind_gt_0.30502884.imiss", header = T)

p7 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count (SNP)") +
    expand_limits(x = c(0,100), y = c(0,6e+05))  +
    labs(title = "After 2nd locus filter + 1st indv filter")
p7


p8 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)", y = "count (indv)", title = "") +
    expand_limits(x = c(0,60), y = c(0,30))  +
    scale_x_continuous(breaks = c(0,15,30,45,60)) +
    geom_vline(xintercept = quantile(imiss$F_MISS, c(0.99))*100, color = "red") +
    annotate(geom = "text", x = (quantile(imiss$F_MISS, c(0.99))*100)+1, hjust = 0, y = 20, label = "99th %") 
p8


#############################################
#Fourth pass -- after remove top 2 indv (99th percentile)


lmiss = read.table("data/Nf/filtering_intermediates/NA_filter.loc_gt_0.3606560.ind_gt_0.23454579.lmiss", header = T)
imiss = read.table("data/Nf/filtering_intermediates/NA_filter.loc_gt_0.3606560.ind_gt_0.23454579.imiss", header = T)

head(lmiss)
head(imiss)

p9 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count (SNP)") +
    geom_vline(xintercept = quantile(lmiss$F_MISS, c(0.95))*100, color = "red") +
    annotate(geom = "text", x = (quantile(lmiss$F_MISS, c(0.95))*100)+1, hjust = 0, y = 3e+05, label = "95th %") +
    expand_limits(x = c(0,100), y = c(0,6e+05))  +
    labs(title = "After 2nd locus filter + 2nd indv filter")
p9


p10 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)", y = "count (indv)", title = "") +
    scale_x_continuous(breaks = c(0,15,30,45,60)) +
    expand_limits(x = c(0,60), y = c(0,30)) 
p10

#############################################
#Fifth pass -- after locus (95th percentile)


lmiss = read.table("data/Nf/filtering_intermediates/NA_filter.loc_gt_0.2750000.ind_gt_0.23454579.lmiss", header = T)
imiss = read.table("data/Nf/filtering_intermediates/NA_filter.loc_gt_0.2750000.ind_gt_0.23454579.imiss", header = T)

p11 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count (SNP)") +
    expand_limits(x = c(0,100), y = c(0,6e+05))  +
    labs(title = "After 3rd locus filter + 2nd indv filter")
p11


p12 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)", y = "count (indv)", title = "") +
    scale_x_continuous(breaks = c(0,15,30,45,60)) +
    geom_vline(xintercept = quantile(imiss$F_MISS, c(0.99))*100, color = "red") +
    annotate(geom = "text", x = (quantile(imiss$F_MISS, c(0.99))*100)+1, hjust = 0, y = 20, label = "99th %") +
    expand_limits(x = c(0,60), y = c(0,30)) 
p12


#############################################
#Sixth pass -- after indv (99th percentile)


lmiss = read.table("data/Nf/filtering_intermediates/NA_filter.loc_gt_0.2750000.ind_gt_0.17720253.lmiss", header = T)
imiss = read.table("data/Nf/filtering_intermediates/NA_filter.loc_gt_0.2750000.ind_gt_0.17720253.imiss", header = T)

p13 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count (SNP)") +
    expand_limits(x = c(0,100), y = c(0,6e+05))  +
    #geom_vline(xintercept = quantile(lmiss$F_MISS, c(0.95))*100, color = "red") +
    geom_vline(xintercept = 25, color = "red") +
    annotate(geom = "text", x = 25+1, hjust = 0, y = 3e+05, label = ">25%") +
    labs(title = "After 3rd locus filter + 3rd indv filter")
p13


p14 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)", y = "count (indv)", title = "") +
    scale_x_continuous(breaks = c(0,15,30,45,60)) +
    expand_limits(x = c(0,60), y = c(0,30)) 
p14


#############################################
#Seventh pass -- after locus (> 25%)


lmiss = read.table("data/Nf/filtering_intermediates/NA_filter.iterative_plus_hard_filter.loc_gt_0.25.ind_gt_0.17720253.lmiss", header = T)
imiss = read.table("data/Nf/filtering_intermediates/NA_filter.iterative_plus_hard_filter.loc_gt_0.25.ind_gt_0.17720253.imiss", header = T)

p15 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count (SNP)") +
    expand_limits(x = c(0,100), y = c(0,6e+05))  +
    labs(title = "After 4th locus filter + 3rd indv filter")
p15

p16 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)", y = "count (indv)", title = "") +
    scale_x_continuous(breaks = c(0,15,30,45,60)) +
    geom_vline(xintercept = 15, color = "red") +
    #annotate(geom = "text", x = (quantile(imiss$F_MISS, c(0.99))*100)+5, y = 20, label = "99th %") +
    annotate(geom = "text", x = 15+1, hjust = 0, y = 20, label = ">15%") +
    expand_limits(x = c(0,60), y = c(0,30)) 
p16

#############################################
#Eighth pass -- after indv (>15%)


lmiss = read.table("data/Nf/filtering_intermediates/NA_filter.iterative_plus_hard_filter.loc_gt_0.25.ind_gt_0.15.lmiss", header = T)
imiss = read.table("data/Nf/filtering_intermediates/NA_filter.iterative_plus_hard_filter.loc_gt_0.25.ind_gt_0.15.imiss", header = T)

p17 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count (SNP)") +
    expand_limits(x = c(0,100), y = c(0,6e+05))  +
    labs(title = "After 4th locus filter + 4th indv filter")
p17

p18 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)", y = "count (indv)", title = "") +
    scale_x_continuous(breaks = c(0,15,30,45,60)) +
    expand_limits(x = c(0,60), y = c(0,30)) 
p18


pdf("figures/quality_filtering/Nf/Nf.lmiss_imiss_distribution.three_pass_sequential_filtering_plus_hard_filter.pdf", height = 16, width = 8)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18, ncol = 2, widths = c(0.5,0.5))
dev.off()


lmiss_iterations_95th = c(0.5242, 0.3606560, 0.2750000, 0.220339 )
plot(lmiss_iterations_95th)
imiss_iterations_99th = c(0.30502884, 0.23454579, 0.17720253, 0.11892679)
plot(imiss_iterations_99th)

lmiss_removed_counts = c(65044, 59085, 56500, 53330)
plot(lmiss_removed_counts)

plot(lmiss_removed_counts~lmiss_iterations_95th)

sum(lmiss_removed_counts)
