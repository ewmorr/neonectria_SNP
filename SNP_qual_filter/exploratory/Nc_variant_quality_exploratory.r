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
    geom_histogram() +
    scale_x_continuous(limits = c(0,7500)) +
    theme_bw()

#QUAL
p5 = ggplot(vcf.qual, aes(x = QUAL)) +
    geom_histogram() +
    #scale_x_continuous(limits = c(0,7500)) +
    theme_bw()
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
) #the ids include both the Nf and Nd samples, so that's fine
second_set.ids = second_set.ids[!second_set.ids %in% third_set.ids]
fourth_set.ids = paste0("NG", 170:196)

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
ind_stats$lib = vector(mode = "character", length = nrow(ind_stats))
ind_stats[row.names(ind_stats) %in% first_set.ids, "lib"] = "lib_1"
ind_stats[row.names(ind_stats) %in% second_set.ids, "lib"] = "lib_2"
ind_stats[row.names(ind_stats) %in% third_set.ids, "lib"] = "lib_3"
ind_stats[row.names(ind_stats) %in% fourth_set.ids, "lib"] = "lib_4"

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

p5 = ggplot(ind_stats, aes(x = lib, y = ind_dp_mean)) +
    geom_boxplot() +
    theme_classic() +
    labs(y = "Individual mean DP", x = "Sequencing library") 
p5

summary(aov(ind_dp_mean ~ lib, data = ind_stats))
# p = 2.42e-08
TukeyHSD(aov(ind_dp_mean ~ lib, data = ind_stats))
#lib 1 dif than libs 2, 3, 4; libs 2 and 3 dif (p = 0.00009), lib 4 dif than 2,3
# i.e., all are different

p6 = ggplot(ind_stats, aes(x = lib, y = ind_na*100)) +
    geom_boxplot() +
    theme_classic() +
    labs(y = "Individual % missing", x = "Sequencing library") 
p6

summary(aov(ind_na ~ lib, data = ind_stats))
# p = 3.69e-06
TukeyHSD(aov(ind_na ~ lib, data = ind_stats))
#lib 1 dif than 2,3,4, no difs among the rest

pdf("figures/quality_filtering/Nc/locus.NA_DP.pdf", width = 10, height = 6)
p1
p3
dev.off()

pdf("figures/quality_filtering/Nc/individual.NA_DP.pdf", width = 16, height = 6)
p2
p4
dev.off()

pdf("figures/quality_filtering/Nc/library.ind_NA_DP.pdf", width = 10, height = 6)
p5
p6
dev.off()

#########################
# there is a not a bimodal distribution in locus mean DP (but there is for Nf)
# it's possible that it's due to library effects.
# split the dp matrix into lib one and lib 2-3 and recalculate locus DP
dp.one = dp[,colnames(dp) %in% first_set.ids]
dp.two = dp[,colnames(dp) %in% second_set.ids]
dp.three = dp[,colnames(dp) %in% third_set.ids]
dp.four = dp[,colnames(dp) %in% fourth_set.ids]

gt.one = gt[,colnames(gt) %in% first_set.ids]
gt.two = gt[,colnames(gt) %in% second_set.ids]
gt.three = gt[,colnames(gt) %in% third_set.ids]
gt.four = gt[,colnames(gt) %in% fourth_set.ids]

dp.one[is.na(gt.one)] = NA
dp.two[is.na(gt.two)] = NA
dp.three[is.na(gt.three)] = NA
dp.four[is.na(gt.four)] = NA


locus_dp_mean.libOne = rowSums(dp.one, na.rm = T)/ncol(dp.one)
locus_dp_mean.libTwo = rowSums(dp.two, na.rm = T)/ncol(dp.two)
locus_dp_mean.libThree = rowSums(dp.three, na.rm = T)/ncol(dp.three)
locus_dp_mean.libFour = rowSums(dp.four, na.rm = T)/ncol(dp.four)

head(locus_dp_mean.libOne)

range(locus_dp_mean.libOne)
range(locus_dp_mean.libTwo)
range(locus_dp_mean.libThree)
range(locus_dp_mean.libFour)

median(locus_dp_mean.libOne)
median(locus_dp_mean.libTwo)
median(locus_dp_mean.libThree)
median(locus_dp_mean.libFour)

mean(locus_dp_mean.libOne)
mean(locus_dp_mean.libTwo)
mean(locus_dp_mean.libThree)
mean(locus_dp_mean.libFour)

sd(locus_dp_mean.libOne)
sd(locus_dp_mean.libTwo)
sd(locus_dp_mean.libThree)
sd(locus_dp_mean.libFour)

mean(locus_dp_mean.libOne)+(2*sd(locus_dp_mean.libOne))
mean(locus_dp_mean.libTwo)+(2*sd(locus_dp_mean.libTwo))
mean(locus_dp_mean.libThree)+(2*sd(locus_dp_mean.libThree))
mean(locus_dp_mean.libFour)+(2*sd(locus_dp_mean.libFour))

mean(locus_dp_mean.libOne)-(1*sd(locus_dp_mean.libOne))
mean(locus_dp_mean.libTwo)-(1*sd(locus_dp_mean.libTwo))
mean(locus_dp_mean.libThree)-(1*sd(locus_dp_mean.libThree))
mean(locus_dp_mean.libFour)-(1*sd(locus_dp_mean.libFour))


quantile(locus_dp_mean.libOne, probs = c(0.05,0.1,0.9,0.95))
quantile(locus_dp_mean.libTwo, probs = c(0.05,0.1,0.9,0.95))
quantile(locus_dp_mean.libThree, probs = c(0.05,0.1,0.9,0.95))
quantile(locus_dp_mean.libFour, probs = c(0.05,0.1,0.9,0.95))

#DP filters
lib_means = c(mean(locus_dp_mean.libOne),
              mean(locus_dp_mean.libTwo),
              mean(locus_dp_mean.libThree),
              mean(locus_dp_mean.libFour)
)

ceiling(lib_means/4)
floor(lib_means*3)

p7 = ggplot(
    data.frame(locus_dp_mean = locus_dp_mean.libOne), 
    aes(x = locus_dp_mean)
) +
    geom_bar() +
    theme_classic() +
    labs(x = "Locus mean DP\n(values of above 100 excluded)", title = "Library one") +
    scale_x_continuous(limits = c(0,100), breaks = c(0,5,10,15,20,25,50,50,75,100))
p7

p8 = ggplot(
    data.frame(locus_dp_mean = locus_dp_mean.libTwo), 
    aes(x = locus_dp_mean)
) +
    geom_bar() +
    theme_classic() +
    labs(x = "Locus mean DP\n(values of above 100 excluded)", title = "Library two") +
    scale_x_continuous(limits = c(0,100), breaks = c(0,5,10,15,20,25,30,35,50,50,75,100))
p8

p9 = ggplot(
    data.frame(locus_dp_mean = locus_dp_mean.libThree), 
    aes(x = locus_dp_mean)
) +
    geom_bar() +
    theme_classic() +
    labs(x = "Locus mean DP\n(values of above 100 excluded)", title = "Library three") +
    scale_x_continuous(limits = c(0,150), breaks = c(0,25,50,50,75,100,125))
p9

p10 = ggplot(
    data.frame(locus_dp_mean = locus_dp_mean.libFour), 
    aes(x = locus_dp_mean)
) +
    geom_bar() +
    theme_classic() +
    labs(x = "Locus mean DP\n(values of above 100 excluded)", title = "Library four") +
    scale_x_continuous(limits = c(0,150), breaks = c(0,25,50,50,75,100,125))
p10

pdf("figures/quality_filtering/Nc/library.locus_DP.pdf", width = 10, height = 6)
p7
p8
p9
p10
dev.off()

#no bimodal distribution in libs
