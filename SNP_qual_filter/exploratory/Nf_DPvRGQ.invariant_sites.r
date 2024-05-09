library(ggplot2)
library(data.table)

#read invariant_sites metrics (geno x indv)
#missingness = fread("~/Nf_invariant_sites_GVCF/NAvREF.invariant_sites.genoXindv.subsample.txt")
#DP = fread("~/Nf_invariant_sites_GVCF/DP.invariant_sites.genoXindv.subsample.txt")
#RGQ = fread("~/Nf_invariant_sites_GVCF/RGQ.invariant_sites.genoXindv.subsample.txt")

#these are already randomly sampled from the full dataset with shuf
missingness = fread("data/Nf_SPANDx_all_seqs/NAvREF.invariant_sites.genoXindv.subsample.txt")
DP = fread("data/Nf_SPANDx_all_seqs/DP.invariant_sites.genoXindv.subsample.txt")
RGQ = fread("data/Nf_SPANDx_all_seqs/RGQ.invariant_sites.genoXindv.subsample.txt")

dat.all = data.frame(
    miss = missingness$V2,
    DP = DP$V1,
    RGQ = RGQ$V1
)

#rm(list = c("missingness", "DP", "RGQ"))
#rand_df <- dat[sample(nrow(dat), size=5*10^5), ]

p1 = ggplot(dat.all, aes(y = DP+1, x = RGQ, color = miss)) +
    geom_point(alpha = 0.15, position = position_jitter(width = 0.4)) +
    scale_color_manual(
        values = c("." = "red", "0" = "black"),
        labels = c("." = "NA", "0" = "REF")
    ) +
    scale_y_log10(breaks = c(1, 10, 100, 500, 1000)) +
    theme_bw()

pdf("figures/DPvRGQ.pdf", width = 8, height = 6)
p1
dev.off()


