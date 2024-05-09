library(vcfR)
library(tidyr)
library(ggplot2)
source("R_scripts/ggplot_theme.txt")

#filtered VCF
vcf <- read.vcfR("data/Nf_SPANDx_all_seqs/out.vcf", verbose = FALSE)
#dna <- ape::read.dna("data/Nf_genome/ref.fasta", format = "fasta")
#gff <- read.table("data/Nf_genome/makerFINAL.all.gff", sep="\t", quote = "")
#rm(dna)
#rm(gff)

#the following pulls in only one tig
#instead let's just try the vcf
#could instead process one tig at a time but we really want dataset wide metrics for filtering here
#chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)


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

pdf("figures/quality_filtering/genotype_qual.SNP.pdf", width = 6, height = 4)
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
third_set.ids = paste0("NG", 170:196)

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
# p < 2e-16
TukeyHSD(aov(ind_dp_mean ~ lib, data = ind_stats))
#lib 1 dif than libs 2 and 3; libs 2 and 3 not dif (p = 0.3)

p6 = ggplot(ind_stats, aes(x = lib, y = ind_na*100)) +
    geom_boxplot() +
    theme_classic() +
    labs(y = "Individual % missing", x = "Sequencing library") 
p6

summary(aov(ind_na ~ lib, data = ind_stats))
# p = 0.688

pdf("figures/quality_filtering/locus.NA_DP.pdf", width = 10, height = 6)
p1
p3
dev.off()

pdf("figures/quality_filtering/individual.NA_DP.pdf", width = 16, height = 6)
p2
p4
dev.off()

pdf("figures/quality_filtering/library.ind_NA_DP.pdf", width = 10, height = 6)
p5
p6
dev.off()

#########################
# there is a bimodal distribution in locus mean DP
# it's possible that it's due to library effects.
# split the dp matrix into lib one and lib 2-3 and recalculate locus DP
dp.one = dp[,colnames(dp) %in% first_set.ids]
dp.two = dp[,colnames(dp) %in% second_set.ids]
dp.three = dp[,colnames(dp) %in% third_set.ids]

locus_dp_mean.libOne = rowSums(dp.one, na.rm = T)/ncol(dp.one)
locus_dp_mean.libTwo = rowSums(dp.two, na.rm = T)/ncol(dp.two)
locus_dp_mean.libThree = rowSums(dp.three, na.rm = T)/ncol(dp.three)

head(locus_dp_mean.libOne)

range(locus_dp_mean.libOne)
range(locus_dp_mean.libTwo)
range(locus_dp_mean.libThree)

median(locus_dp_mean.libOne)
median(locus_dp_mean.libTwo)
median(locus_dp_mean.libThree)

mean(locus_dp_mean.libOne)
mean(locus_dp_mean.libTwo)
mean(locus_dp_mean.libThree)


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
    scale_x_continuous(limits = c(0,100), breaks = c(0,25,50,50,75,100))
p9

pdf("figures/quality_filtering/library.locus_NA.pdf", width = 10, height = 6)
p7
p8
p9
dev.off()














#################################
#OLD stuff
##############


#set dp to NA where genotype is na
dp[is.na(gt)] = NA
sum(dp == 0, na.rm = T)
sum(dp == 1)


dp[dp == 0] = NA
apply(dp, 2, median, na.rm = T)
apply(dp, 1, median, na.rm = T)
median(dp, na.rm = T)
mean(dp, na.rm = T)

dpf <- pivot_longer(data.frame(dp), names_to = "Sample", values_to = "Depth", cols = everything()) #DAMN YOU HADLEY!!!!
dpf <- dpf[ dpf$Depth > 0,]

p1 <- ggplot(dpf, aes(x=Sample, y=Depth)) +
    geom_violin(fill="#C0C0C0", adjust=1.0, scale = "count", trim=TRUE) +
    stat_summary(fun=median, geom="point", shape=23, size=2) +
    my_gg_theme +
    scale_y_continuous(trans=scales::log2_trans(),
        breaks=c(1, 10, 100, 800),
        minor_breaks=c(1:10, 2:10*10, 2:8*100)
    ) +
    labs(x = "samples") +
    theme(
        axis.text.x = element_blank()
    )
p1

pdf("depth_by_samples.filtered.PASS.pdf", width = 12, height = 6)
p1
dev.off()

#median/mean coverage
median_cov = apply(dp, 2, median, na.rm = T)
mean_cov = apply(dp, 2, mean, na.rm = T)

mean(mean_cov)

cov_mat = matrix(c("median", median(median_cov), "mean", mean(mean_cov) ), nrow = 2)
write.table(data.frame(cov_mat), "VCF.filtered.PASS.coverage.txt", quote = F, sep = "\t", row.names = F, col.names = F)

#unfiltered VCF
vcf <- read.vcfR("~/Nf_SPANDx_all_seqs/Outputs/Master_vcf/out.vcf", verbose = FALSE)
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)

dpf <- pivot_longer(data.frame(dp), names_to = "Sample", values_to = "Depth", cols = everything()) #DAMN YOU HADLEY!!!!
dpf <- dpf[ dpf$Depth > 0,]

p1 <- ggplot(dpf, aes(x=Sample, y=Depth)) +
geom_violin(fill="#C0C0C0", adjust=1.0, scale = "count", trim=TRUE) +
stat_summary(fun=median, geom="point", shape=23, size=2) +
my_gg_theme +
scale_y_continuous(trans=scales::log2_trans(),
breaks=c(1, 10, 100, 800),
minor_breaks=c(1:10, 2:10*10, 2:8*100)
) +
labs(x = "samples") +
theme(
axis.text.x = element_blank()
)

pdf("depth_by_samples.unfiltered.pdf", width = 12, height = 6)
p1
dev.off()

#median/mean coverage
median_cov = apply(dp, 2, median, na.rm = T)
mean_cov = apply(dp, 2, mean, na.rm = T)

mean(mean_cov)

cov_mat = matrix(c("median", median(median_cov), "mean", mean(mean_cov) ), nrow = 2)
write.table(data.frame(cov_mat), "VCF.unfiltered.coverage.txt", quote = F, sep = "\t", row.names = F, col.names = F)
