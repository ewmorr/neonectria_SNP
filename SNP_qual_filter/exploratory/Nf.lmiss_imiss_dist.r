library(ggplot2)
library(gridExtra)

lmiss = read.table("data/Nf/all_libs.DPGQ_filter.lmiss", header = T)
imiss = read.table("data/Nf/all_libs.DPGQ_filter.imiss", header = T)

head(lmiss)
head(imiss)

p1 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count") +
    expand_limits(x = c(0,100), y = c(0,5e+05)) 
p1
# looks approximately half-normal
# SD of half-normal as described here https://stats.stackexchange.com/questions/535043/standard-deviation-for-a-half-normal-distribution

hn_sd = function(x){
    return( sqrt( sum(x^2)/length(x) ) )
}

lmiss_sd = hn_sd(lmiss$F_MISS)
2*lmiss_sd
quantile(lmiss$F_MISS, c(0.25, 0.5, 0.75, 0.95))
quantile(lmiss$F_MISS, c(0.80, 0.85, 0.90, 0.95))

p1 = p1 + #geom_vline(xintercept = lmiss_sd*2*100, color = "blue") +
    geom_vline(xintercept = quantile(lmiss$F_MISS, c(0.95))*100, color = "red") +
    #annotate(geom = "text", x = (lmiss_sd*2*100)-7, y = 3e+05, label = "2*SD") +
    annotate(geom = "text", x = (quantile(lmiss$F_MISS, c(0.95))*100)+7, y = 3e+05, label = "95th %") +
    labs(title = "no missingness filter")

p2 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)") +
    expand_limits(x = c(0,60), y = c(0,27))  +
    labs(title = "")
p2

imiss_m = mean(imiss$F_MISS)
imiss_sd = sd(imiss$F_MISS)
imiss_m
imiss_sd
imiss_m+(2*imiss_sd)

quantile(imiss$F_MISS, c(0.25, 0.5, 0.75, 0.95))
quantile(imiss$F_MISS, c(0.25, 0.5, 0.75, 0.99))

#p2 = p2 + geom_vline(xintercept = (imiss_m+(2*imiss_sd))*100, color = "blue") +
#    geom_vline(xintercept = quantile(imiss$F_MISS, c(0.95))*100, color = "red") +
#    labs(title = "")


pdf("figures/quality_filtering/Nf/Nf.lmiss_miss.dist.pdf", height = 5, width = 10)
grid.arrange(p1,p2,ncol = 2)
dev.off()

###########################################
#Import new tables and rerun above

lmiss = read.table("data/Nf/NA_filter.loc_gt_0.524194.lmiss", header = T)
imiss = read.table("data/Nf/NA_filter.loc_gt_0.524194.imiss", header = T)

head(lmiss)
head(imiss)

p1 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count") +
    expand_limits(x = c(0,100), y = c(0,5e+05)) 
p1
# looks approximately half-normal
# SD of half-normal as described here https://stats.stackexchange.com/questions/535043/standard-deviation-for-a-half-normal-distribution

hn_sd = function(x){
    return( sqrt( sum(x^2)/length(x) ) )
}

lmiss_sd = hn_sd(lmiss$F_MISS)
2*lmiss_sd
quantile(lmiss$F_MISS, c(0.25, 0.5, 0.75, 0.95))
quantile(lmiss$F_MISS, c(0.80, 0.85, 0.90, 0.95))

p1 = p1 + geom_vline(xintercept = lmiss_sd*2*100, color = "blue") +
    geom_vline(xintercept = quantile(lmiss$F_MISS, c(0.95))*100, color = "red") +
    annotate(geom = "text", x = (lmiss_sd*2*100)-7, y = 3e+05, label = "2*SD") +
    annotate(geom = "text", x = (quantile(lmiss$F_MISS, c(0.95))*100)+7, y = 3e+05, label = "95th %") +
    labs(title = "1st locus missingness filter")

p2 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)") +
    expand_limits(x = c(0,60), y = c(0,27)) 
p2

imiss_m = mean(imiss$F_MISS)
imiss_sd = sd(imiss$F_MISS)
imiss_m
imiss_sd
imiss_m+(2*imiss_sd)

quantile(imiss$F_MISS, c(0.25, 0.5, 0.75, 0.95))
quantile(imiss$F_MISS, c(0.25, 0.5, 0.75, 0.99))

p2 = p2 + geom_vline(xintercept = (imiss_m+(2*imiss_sd))*100, color = "blue") +
    geom_vline(xintercept = quantile(imiss$F_MISS, c(0.95))*100, color = "red") +
    labs(title = "")


pdf("figures/quality_filtering/Nf/Nf.lmiss_miss.dist.1st_lNA_filter.pdf", height = 5, width = 10)
grid.arrange(p1,p2,ncol = 2)
dev.off()

#############################################
#Second pass -- removed top indv NA (NG 139)


lmiss = read.table("data/Nf/NA_filter.loc_gt_0.524194.ind_gt_0.30502884.lmiss", header = T)
imiss = read.table("data/Nf/NA_filter.loc_gt_0.524194.ind_gt_0.30502884.imiss", header = T)

head(lmiss)
head(imiss)

p1 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count") +
    expand_limits(x = c(0,100), y = c(0,5e+05)) 
p1

lmiss_sd = hn_sd(lmiss$F_MISS)
2*lmiss_sd
quantile(lmiss$F_MISS, c(0.25, 0.5, 0.75, 0.95))
quantile(lmiss$F_MISS, c(0.80, 0.85, 0.90, 0.95))

p1 = p1 + geom_vline(xintercept = lmiss_sd*2*100, color = "blue") +
    geom_vline(xintercept = quantile(lmiss$F_MISS, c(0.95))*100, color = "red") +
    annotate(geom = "text", x = (lmiss_sd*2*100)-7, y = 2e+05, label = "2*SD") +
    annotate(geom = "text", x = (quantile(lmiss$F_MISS, c(0.95))*100)+7, y = 2e+05, label = "95th %") +
    labs(title = "1st locus + 1st indv missingness filter")

p2 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)", y = "count") +
    expand_limits(x = c(0,60), y = c(0,27)) 
p2

imiss_m = mean(imiss$F_MISS)
imiss_sd = sd(imiss$F_MISS)
imiss_m
imiss_sd
imiss_m+(2*imiss_sd)

quantile(imiss$F_MISS, c(0.25, 0.5, 0.75, 0.95))
quantile(imiss$F_MISS, c(0.25, 0.5, 0.75, 0.99))

p2 = p2 + geom_vline(xintercept = (imiss_m+(2*imiss_sd))*100, color = "blue") +
    geom_vline(xintercept = quantile(imiss$F_MISS, c(0.95))*100, color = "red") +
    geom_vline(xintercept = quantile(imiss$F_MISS, c(0.99))*100, color = "green") +
    annotate(geom = "text", x = (quantile(imiss$F_MISS, c(0.99))*100)+2, y = 10, label = "99th %") +
    labs(title = "")


pdf("figures/quality_filtering/Nf/Nf.lmiss_miss.dist.1st_lNA_1st_iNA_filter.pdf", height = 5, width = 10)
grid.arrange(p1,p2,ncol = 2)
dev.off()

#############################################
#Third pass -- remove top 5% loci


lmiss = read.table("data/Nf/NA_filter.loc_gt_0.3606560.ind_gt_0.30502884.lmiss", header = T)
imiss = read.table("data/Nf/NA_filter.loc_gt_0.3606560.ind_gt_0.30502884.imiss", header = T)

head(lmiss)
head(imiss)

p1 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count") +
    expand_limits(x = c(0,100), y = c(0,5e+05)) 
p1

lmiss_sd = hn_sd(lmiss$F_MISS)
2*lmiss_sd
quantile(lmiss$F_MISS, c(0.25, 0.5, 0.75, 0.95))
quantile(lmiss$F_MISS, c(0.80, 0.85, 0.90, 0.95))

p1 = p1 + geom_vline(xintercept = lmiss_sd*2*100, color = "blue") +
    geom_vline(xintercept = quantile(lmiss$F_MISS, c(0.95))*100, color = "red") +
    annotate(geom = "text", x = (lmiss_sd*2*100)-3, y = 2e+05, label = "2*SD") +
    annotate(geom = "text", x = (quantile(lmiss$F_MISS, c(0.95))*100)+3, y = 2e+05, label = "95th %") +
    labs(title = "2nd locus + 1st indv missingness filter")

p2 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)") +
    expand_limits(x = c(0,60), y = c(0,27)) 
p2

imiss_m = mean(imiss$F_MISS)
imiss_sd = sd(imiss$F_MISS)
imiss_m
imiss_sd
imiss_m+(2*imiss_sd)

quantile(imiss$F_MISS, c(0.25, 0.5, 0.75, 0.95))
quantile(imiss$F_MISS, c(0.25, 0.5, 0.75, 0.99))

p2 = p2 + geom_vline(xintercept = (imiss_m+(2*imiss_sd))*100, color = "blue") +
    geom_vline(xintercept = quantile(imiss$F_MISS, c(0.95))*100, color = "red") +
    geom_vline(xintercept = quantile(imiss$F_MISS, c(0.99))*100, color = "green") +
    annotate(geom = "text", x = (quantile(imiss$F_MISS, c(0.99))*100)+2, y = 10, label = "99th %") +
    labs(title = "")
p2

pdf("figures/quality_filtering/Nf/Nf.lmiss_miss.dist.2nd_lNA_1st_iNA_filter.pdf", height = 5, width = 10)
grid.arrange(p1,p2,ncol = 2)
dev.off()

#############################################
#Fourth pass -- remove top 2 indv (99th percentile)


lmiss = read.table("data/Nf/NA_filter.loc_gt_0.3606560.ind_gt_0.23454579.lmiss", header = T)
imiss = read.table("data/Nf/NA_filter.loc_gt_0.3606560.ind_gt_0.23454579.imiss", header = T)

head(lmiss)
head(imiss)

p1 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count") +
    expand_limits(x = c(0,100), y = c(0,5e+05)) 
p1

lmiss_sd = hn_sd(lmiss$F_MISS)
2*lmiss_sd
quantile(lmiss$F_MISS, c(0.25, 0.5, 0.75, 0.95))
quantile(lmiss$F_MISS, c(0.80, 0.85, 0.90, 0.95))

p1 = p1 + geom_vline(xintercept = lmiss_sd*2*100, color = "blue") +
    geom_vline(xintercept = quantile(lmiss$F_MISS, c(0.95))*100, color = "red") +
    annotate(geom = "text", x = (lmiss_sd*2*100)-3, y = 2e+05, label = "2*SD") +
    annotate(geom = "text", x = (quantile(lmiss$F_MISS, c(0.95))*100)+3, y = 2e+05, label = "95th %") +
    labs(title = "2nd locus + 2nd indv missingness filter")

p2 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)", y = "count") +
    expand_limits(x = c(0,60), y = c(0,27)) 
p2

imiss_m = mean(imiss$F_MISS)
imiss_sd = sd(imiss$F_MISS)
imiss_m
imiss_sd
imiss_m+(2*imiss_sd)

quantile(imiss$F_MISS, c(0.25, 0.5, 0.75, 0.95))
quantile(imiss$F_MISS, c(0.25, 0.5, 0.75, 0.99))

p2 = p2 + geom_vline(xintercept = (imiss_m+(2*imiss_sd))*100, color = "blue") +
    geom_vline(xintercept = quantile(imiss$F_MISS, c(0.95))*100, color = "red") +
    geom_vline(xintercept = quantile(imiss$F_MISS, c(0.99))*100, color = "green") +
    annotate(geom = "text", x = (quantile(imiss$F_MISS, c(0.99))*100)+2, y = 10, label = "99th %") +
    labs(title = "")
p2

pdf("figures/quality_filtering/Nf/Nf.lmiss_miss.dist.2nd_lNA_2nd_iNA_filter.pdf", height = 5, width = 10)
grid.arrange(p1,p2,ncol = 2)
dev.off()

#############################################
#Fifth pass -- locus (95th percentile)


lmiss = read.table("data/Nf/NA_filter.loc_gt_0.2750000.ind_gt_0.23454579.lmiss", header = T)
imiss = read.table("data/Nf/NA_filter.loc_gt_0.2750000.ind_gt_0.23454579.imiss", header = T)

head(lmiss)
head(imiss)

p1 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count") +
    expand_limits(x = c(0,100), y = c(0,5e+05)) 
p1

lmiss_sd = hn_sd(lmiss$F_MISS)
2*lmiss_sd
quantile(lmiss$F_MISS, c(0.25, 0.5, 0.75, 0.95))
quantile(lmiss$F_MISS, c(0.80, 0.85, 0.90, 0.95))

p1 = p1 + geom_vline(xintercept = lmiss_sd*2*100, color = "blue") +
    geom_vline(xintercept = quantile(lmiss$F_MISS, c(0.95))*100, color = "red") +
    annotate(geom = "text", x = (lmiss_sd*2*100)-3, y = 2e+05, label = "2*SD") +
    annotate(geom = "text", x = (quantile(lmiss$F_MISS, c(0.95))*100)+3, y = 2e+05, label = "95th %") +
    labs(title = "3rd locus + 2nd indv missingness filter")

p2 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)", y = "count") +
    expand_limits(x = c(0,60), y = c(0,27)) 
p2

imiss_m = mean(imiss$F_MISS)
imiss_sd = sd(imiss$F_MISS)
imiss_m
imiss_sd
imiss_m+(2*imiss_sd)

quantile(imiss$F_MISS, c(0.25, 0.5, 0.75, 0.95))
quantile(imiss$F_MISS, c(0.25, 0.5, 0.75, 0.99))

p2 = p2 + geom_vline(xintercept = (imiss_m+(2*imiss_sd))*100, color = "blue") +
    geom_vline(xintercept = quantile(imiss$F_MISS, c(0.95))*100, color = "red") +
    geom_vline(xintercept = quantile(imiss$F_MISS, c(0.99))*100, color = "green") +
    annotate(geom = "text", x = (quantile(imiss$F_MISS, c(0.99))*100)+2, y = 10, label = "99th %") +
    labs(title = "")
p2

pdf("figures/quality_filtering/Nf/Nf.lmiss_miss.dist.3rd_lNA_2nd_iNA_filter.pdf", height = 5, width = 10)
grid.arrange(p1,p2,ncol = 2)
dev.off()



#############################################
#Sixth pass -- indv (99th percentile)


lmiss = read.table("data/Nf/NA_filter.loc_gt_0.2750000.ind_gt_0.17720253.lmiss", header = T)
imiss = read.table("data/Nf/NA_filter.loc_gt_0.2750000.ind_gt_0.17720253.imiss", header = T)

head(lmiss)
head(imiss)

p1 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count") +
    expand_limits(x = c(0,100), y = c(0,5e+05)) 
p1

lmiss_sd = hn_sd(lmiss$F_MISS)
2*lmiss_sd
quantile(lmiss$F_MISS, c(0.25, 0.5, 0.75, 0.95))
quantile(lmiss$F_MISS, c(0.80, 0.85, 0.90, 0.95))

p1 = p1 + geom_vline(xintercept = lmiss_sd*2*100, color = "blue") +
    geom_vline(xintercept = quantile(lmiss$F_MISS, c(0.95))*100, color = "red") +
    annotate(geom = "text", x = (lmiss_sd*2*100)-3, y = 2e+05, label = "2*SD") +
    annotate(geom = "text", x = (quantile(lmiss$F_MISS, c(0.95))*100)+3, y = 2e+05, label = "95th %") +
    labs(title = "3rd locus + 3rd indv missingness filter")

p2 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)", y = "count") +
    expand_limits(x = c(0,60), y = c(0,27)) 
p2

imiss_m = mean(imiss$F_MISS)
imiss_sd = sd(imiss$F_MISS)
imiss_m
imiss_sd
imiss_m+(2*imiss_sd)

quantile(imiss$F_MISS, c(0.25, 0.5, 0.75, 0.95))
quantile(imiss$F_MISS, c(0.25, 0.5, 0.75, 0.99))

p2 = p2 + geom_vline(xintercept = (imiss_m+(2*imiss_sd))*100, color = "blue") +
    geom_vline(xintercept = quantile(imiss$F_MISS, c(0.95))*100, color = "red") +
    geom_vline(xintercept = quantile(imiss$F_MISS, c(0.99))*100, color = "green") +
    annotate(geom = "text", x = (quantile(imiss$F_MISS, c(0.99))*100)+2, y = 10, label = "99th %") +
    labs(title = "")
p2

pdf("figures/quality_filtering/Nf/Nf.lmiss_miss.dist.3rd_lNA_3rd_iNA_filter.pdf", height = 5, width = 10)
grid.arrange(p1,p2,ncol = 2)
dev.off()


lmiss_iterations_95th = c(0.5242, 0.36290300, 0.3606560, 0.2786890, 0.2750000, 0.225000, 0.220339 )
plot(lmiss_iterations_95th)
imiss_iterations_99th = c(0.32507473, 0.30502884, 0.24944936, 0.23454579, 0.18945936 , 0.17720253, 0.12870211)
plot(imiss_iterations_99th)

#############################################
#Seventh pass -- locus (95th percentile)


lmiss = read.table("data/Nf/NA_filter.loc_gt_0.220339.ind_gt_0.17720253.lmiss", header = T)
imiss = read.table("data/Nf/NA_filter.loc_gt_0.220339.ind_gt_0.17720253.imiss", header = T)

p1 = ggplot(lmiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing individuals per locus (%)", y = "count") +
    expand_limits(x = c(0,100), y = c(0,5e+05)) 
p1

lmiss_sd = hn_sd(lmiss$F_MISS)
2*lmiss_sd
quantile(lmiss$F_MISS, c(0.25, 0.5, 0.75, 0.95))
quantile(lmiss$F_MISS, c(0.80, 0.85, 0.90, 0.95))

p1 = p1 + geom_vline(xintercept = lmiss_sd*2*100, color = "blue") +
    geom_vline(xintercept = quantile(lmiss$F_MISS, c(0.95))*100, color = "red") +
    annotate(geom = "text", x = (lmiss_sd*2*100)-3, y = 2e+05, label = "2*SD") +
    annotate(geom = "text", x = (quantile(lmiss$F_MISS, c(0.95))*100)+3, y = 2e+05, label = "95th %") +
    labs(title = "4th locus + 3rd indv missingness filter")

p2 = ggplot(imiss, aes(x = F_MISS*100)) +
    geom_histogram(bins = 50) +
    theme_classic() +
    labs(x = "Missing loci per indidual (%)", y = "count") +
    expand_limits(x = c(0,60), y = c(0,27)) 
p2

imiss_m = mean(imiss$F_MISS)
imiss_sd = sd(imiss$F_MISS)
imiss_m
imiss_sd
imiss_m+(2*imiss_sd)

quantile(imiss$F_MISS, c(0.25, 0.5, 0.75, 0.95))
quantile(imiss$F_MISS, c(0.25, 0.5, 0.75, 0.99))

p2 = p2 + geom_vline(xintercept = (imiss_m+(2*imiss_sd))*100, color = "blue") +
    geom_vline(xintercept = quantile(imiss$F_MISS, c(0.95))*100, color = "red") +
    geom_vline(xintercept = quantile(imiss$F_MISS, c(0.99))*100, color = "green") +
    annotate(geom = "text", x = (quantile(imiss$F_MISS, c(0.99))*100)+2, y = 10, label = "99th %") +
    labs(title = "")
p2

pdf("figures/quality_filtering/Nf/Nf.lmiss_miss.dist.4th_lNA_3rd_iNA_filter.pdf", height = 5, width = 10)
grid.arrange(p1,p2,ncol = 2)
dev.off()


lmiss_iterations_95th = c(0.5242, 0.36290300, 0.3606560, 0.2786890, 0.2750000, 0.225000, 0.220339 )
plot(lmiss_iterations_95th)
imiss_iterations_99th = c(0.32507473, 0.30502884, 0.24944936, 0.23454579, 0.18945936 , 0.17720253, 0.12870211)
plot(imiss_iterations_99th)

#############################################
#Eighth pass -- locus (95th percentile)


lmiss = read.table("data/Nf/NA_filter.loc_gt_0.220339.ind_gt_0.11892679.lmiss", header = T)
imiss = read.table("data/Nf/NA_filter.loc_gt_0.220339.ind_gt_0.11892679.imiss", header = T)


lmiss_sd = hn_sd(lmiss$F_MISS)
2*lmiss_sd
quantile(lmiss$F_MISS, c(0.25, 0.5, 0.75, 0.95))
quantile(lmiss$F_MISS, c(0.80, 0.85, 0.90, 0.95))

max(lmiss$F_MISS)

imiss_m = mean(imiss$F_MISS)
imiss_sd = sd(imiss$F_MISS)
imiss_m
imiss_sd
imiss_m+(2*imiss_sd)

quantile(imiss$F_MISS, c(0.25, 0.5, 0.75, 0.95))
quantile(imiss$F_MISS, c(0.25, 0.5, 0.75, 0.99))

max(imiss$F_MISS)
