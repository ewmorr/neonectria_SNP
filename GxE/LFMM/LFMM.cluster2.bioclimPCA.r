library(lfmm)
library(dplyr)
library(ggplot2)
library(scales)
#This function for plotting reversed log10 of p vals
reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv,
    log_breaks(base = base),
    domain = c(1e-100, Inf))
}
source("library/ggplot_theme.txt")

#read site metadata
site.info = read.csv("data/sample_metadata/site_info.csv")
site.bioclim = read.csv("data/sample_metadata/Nf.site_bioclim.csv")
site.bioclim = site.bioclim %>% filter(state != "WV")
site.bioclim.PC = read.csv("data/sample_metadata/Nf.sites_bioclim_PC.no_WV.csv", header = T)
site.bioclim
site.bioclim.PC

site_metadata = left_join(
    site.bioclim,
    site.info %>% select(state, lat, lon, duration_infection),
    by = c("state", "lat", "lon")
) %>% unique() %>% left_join(
    ., site.bioclim.PC %>% select(-bio2),
    by = "state"
)



#sample metadata. Needs to be sorted by PED .sampleIDs info
sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
rownames(sample_metadata.Nf) = sample_metadata.Nf$Sequence_label
#PED sample IDs
sample_ids = read.table("data/Nf/final_tables/rm_dups/FINAL_snp.gwas_analyses.cluster2.sampleIDs", header = F)
sample_metadata.Nf.sorted = sample_metadata.Nf[sample_ids$V1,]
nrow(sample_metadata.Nf.sorted)
sample_metadata.site_info = left_join(sample_metadata.Nf.sorted, site_metadata)
is.na(sample_metadata.site_info$bio1) %>% sum()

row_ids = which(sample_ids$V1 %in% sample_metadata.site_info$Sequence_label)

#The genotype data can simply be read in as a matrix (according the docs)
#OR can try loading LEA and using readLfmm()
Y = as.matrix(read.table("data/Nf/final_tables/rm_dups/FINAL_snp.gwas_analyses.cluster2.lfmm", header = F))
Y.filtered = Y[row_ids,]
nrow(Y.filtered)
ncol(Y.filtered)

SNP_pos = read.table("data/Nf/final_tables/rm_dups/FINAL_snp.gwas_analyses.cluster2.recode.map")
nrow(SNP_pos)
SNP_pos = SNP_pos[c(1,4)]
colnames(SNP_pos) = c("scaffold", "position")


###################
###################
###################
#skip this routine
ref_sum = colSums(Y.filtered == 0)
alt_sum = colSums(Y.filtered == 1)

#filter for MAC ge 3
nrow(Y.filtered)
minMAC = 3

sum(ref_sum < minMAC)
# 101504
sum(alt_sum < minMAC)
ncol(Y.filtered)
# 424811
ncol(Y.filtered) - sum(ref_sum < minMAC)
# 323307

which(ref_sum < minMAC) 
which(alt_sum < minMAC)# there are none
rm_cols = which(ref_sum < minMAC)

length(rm_cols)\#101504
ncol(Y.filtered)
length(rm_cols)/ncol(Y.filtered)
# 0.2389392
ncol(Y.filtered)-length(rm_cols)
# 323307
Y.filteredMAC = Y.filtered[,-rm_cols]
SNP_pos.filteredMAC = SNP_pos[-rm_cols,]
ncol(Y.filteredMAC)
nrow(SNP_pos.filteredMAC)

#filter for greater than 25% missing data
#sum(colSums(Y.filteredMAC == 9)/nrow(Y.filteredMAC) > 0.25)
#rm_cols = which(colSums(Y.filteredMAC == 9)/nrow(Y.filteredMAC) > 0.25)
#Y.filteredNA = Y.filteredMAC[,-rm_cols]

plot(colSums(Y == 9)/nrow(Y))
plot(colSums(Y == 0))
plot(colSums(Y == 1))
###################
###################
###################

Y.filtered = Y.filteredMAC
SNP_pos = SNP_pos.filteredMAC



#
#principal components analysis for K
pc = prcomp(Y.filtered)
plot((pc$sdev^2)/sum(pc$sdev^2), xlab = 'PC', ylab = "% variance explained")
(pc$sdev^2)/sum(pc$sdev^2)
points(4,pc$sdev[7]^2, type = "h", lwd = 3, col = "blue")
str(pc)
plot(pc$x[,2] ~ pc$x[,1])
plot(pc$x[,3] ~ pc$x[,1])
plot(pc$x[,3] ~ pc$x[,2])
plot(pc$x[,4] ~ pc$x[,1])

#Y = Y.filteredNA
#######################
#PC1

#variable for test
X = (sample_metadata.site_info$temp.PC1)
#X = (sample_metadata.site_info[,39:43])


class(X)
class(Y.filtered)
sd(X)
mean(X)
#LFMM ridge

mod.lfmm = lfmm_ridge(Y = Y.filtered, X = X, K = 2, lambda = 1) #using K = 2 based on PCA 
str(mod.lfmm)

pv <- lfmm_test(Y = Y.filtered,
X = X,
lfmm = mod.lfmm,
calibrate = "gif")
str(pv)

#Example plots
plot(-log10(pv$calibrated.pvalue),
pch = 19,
cex = .3,
xlab = "SNP", ylab = "-Log P",
col = "grey")

plot((pv$B),
pch = 19,
cex = .3,
xlab = "SNP", ylab = "Effect size",
col = "grey")

#Computing genomic inflation factor (GIF) based on calibrated z-scores (http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html) and Francois et al. 2016
lambda = median(pv$score^2, na.rm = T)/0.456
lambda #1.180159
adj.p.values = pchisq(pv$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
#Note that these calibrated scores are similar as pv$calibrated.pvalue
hist(pv$calibrated.pvalue) #this looks basically perfect
hist(pv$pvalue)
#Try higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv$score^2/1.25, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv$score^2/1.05, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv$score^2/1.05, df = 1, lower = FALSE)
hist(adj.p.values)

#THIS IS SHOWING THAT THE GIF CALIBRATION IN THE ALGORITHM IS MORE CONSERVATIVE THAN LOWER VALUES OF LAMBDA
#However, the lower values have a correct distribution under null model
#Try at GIF = 1.255808
adj.p.values = pchisq(pv$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
length(adj.p.values)

#Read scaffold lengths
scf_lens = read.table("data/Nf/final_tables/rm_dups/scaffold_lengths.csv", sep = ",", header = F)
colnames(scf_lens) = c("scaffold", "length")

#Join with actual positiion and chromosome
nrow(SNP_pos)
ncol(Y.filtered)
pv.with_pos = data.frame(calibrated.p = pv$calibrated.pvalue, effect_size = pv$B, SNP_pos, man.adj.p = adj.p.values) %>% left_join(., scf_lens)

#############################
#Outlier tests
# identify the 95% percentile
my_threshold <- quantile((pv.with_pos )$calibrated.p, 0.025, na.rm = T) #removed the filter by 100000, can do this later for plotting but it does not seem good to do before outlier ID; %>% filter(length > 100000)
# make an outlier column in the data.frame
pv.with_pos <- pv.with_pos %>% mutate(outlier = ifelse(calibrated.p < my_threshold, "outlier", "background"))
#Number of outliers
#note that once we filtered for MAC we have no more NA p values
pv.with_pos %>% group_by(outlier) %>% tally()
pv.with_pos %>% filter(is.na(outlier)) %>% head
#8083
#FDR correction
#This is based on the auto calibartion
pv.with_pos$FDR.p = p.adjust(pv.with_pos$calibrated.p, method = "fdr", n = length(pv.with_pos$calibrated.p))
pv.with_pos <- pv.with_pos %>% mutate(FDR.sig = ifelse(FDR.p < 0.05, "sig", "background"))
pv.with_pos %>% group_by(FDR.sig) %>% tally()

#1092  SNPs identified as significant after FDR correction


#FDR correction
#This is based on the manual GIF adjustment
pv.with_pos$FDR.p.man = p.adjust(pv.with_pos$man.adj.p, method = "fdr", n = length(pv.with_pos$man.adj.p))
pv.with_pos <- pv.with_pos %>% mutate(FDR.sig.man = ifelse(FDR.p.man < 0.05, "sig", "background"))
pv.with_pos %>% group_by(FDR.sig.man) %>% tally()


pv.PC1.with_pos = pv.with_pos
#1094 SNPs identified as significant after FDR correction


pv.PC1.with_pos %>% head()
pv.PC1.with_pos %>% filter(outlier == "outlier") %>% nrow()
#8083 
# but this is not the same thing as sig these are 2.5% outliers

####################
#ggplots

#Basic plot
ggplot(
    pv.PC1.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = calibrated.p)
) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
#scale_x_continuous(labels = fancy_scientific, breaks = c(1, seq(from = 10^6, to = 6*10^6, by = 10^6)) ) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
my_gg_theme +
labs(x = "Position (Mbp)", y = "P value") +
theme(
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 8)
    #axis.text.x = element_text(angle = 85, size = 10, hjust = 1)
)

#Colored by outliers
ggplot(
    pv.PC1.with_pos %>% filter(length > 100000 ), 
    aes(x = position/10^6, y = calibrated.p, color = outlier)
) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black")) +
my_gg_theme +
labs(x = "Position (Mbp)", y = "P value") +
theme(
    strip.text.x = element_blank(),
    axis.text.x = element_text(size = 8),
    legend.position = "bottom"
)

#FDR with auto corrected P (algorithm GIF)
p1 = ggplot(
    pv.PC1.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = calibrated.p, color = FDR.sig)
) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "") +
theme(
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)
p1

p2 = ggplot(
    pv.PC1.with_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man)
) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "") +
theme(
    strip.text.x = element_blank(),
    #axis.text.x = element_text(size = 8)
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)
p2

png("figures/GxE/LFMM/bioclim.temp.PC1.png", width = 1080, height = 240)
p1
dev.off()

#######################
#PC2
#######


#variable for test
X = sample_metadata.site_info$precip.PC1

#LFMM ridge
#mod.lfmm = lfmm_ridge(Y = Y, X = X, K = 2)
mod.lfmm.pc2 = lfmm_ridge(Y = Y.filtered, X = X, K = 2)

pv.pc2 <- lfmm_test(Y = Y.filtered,
X = X,
lfmm = mod.lfmm.pc2,
calibrate = "gif")

#Example plots
plot(-log10(pv.pc2$calibrated.pvalue),
pch = 19,
cex = .3,
xlab = "Probe", ylab = "-Log P",
col = "grey")

plot((pv.pc2$B),
pch = 19,
cex = .3,
xlab = "Probe", ylab = "Effect size",
col = "grey")

#Computing genomic inflation factor (GIF) based on calibrated z-scores (http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html) and Francois et al. 2016
lambda = median(pv.pc2$score^2, na.rm = T)/0.456
lambda #1.294011
adj.p.values = pchisq(pv.pc2$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
#Note that these calibrated scores are similar as pv$calibrated.pvalue
hist(pv.pc2$calibrated.pvalue)
hist(pv.pc2$pvalue)
#IN THIS CASE THE CALCULATED VALUES AREADY LOOK GOOD 

#Try higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.pc2$score^2/1.35, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.pc2$score^2/1.15, df = 1, lower = FALSE)
hist(adj.p.values)

#THIS IS SHOWING THAT THE GIF CALIBRATION IN THE ALGORITHM IS Good
#However, the lower values have a correct distribution under null model
#Try at GIF = lambda
adj.p.values = pchisq(pv.pc2$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)

nrow(SNP_pos)
nrow(scf_lens)
#Join with actual positiion and chromosome
pv.pc2.with_pos = data.frame(calibrated.p = pv.pc2$calibrated.pvalue, effect_size = pv.pc2$B, SNP_pos, man.adj.p = adj.p.values) %>% left_join(., scf_lens)

#############################
#Outlier tests
# identify the 95% percentile
my_threshold <- quantile((pv.pc2.with_pos)$calibrated.p, 0.025, na.rm = T)
# make an outlier column in the data.frame
pv.pc2.with_pos <- pv.pc2.with_pos %>% mutate(outlier = ifelse(calibrated.p < my_threshold, "outlier", "background"))
#Number of outliers
pv.pc2.with_pos %>% group_by(outlier) %>% tally()

#8083

#FDR correction
pv.pc2.with_pos = pv.pc2.with_pos %>% filter(!is.na(calibrated.p))

pv.pc2.with_pos$FDR.p = p.adjust(pv.pc2.with_pos$calibrated.p, method = "fdr", n = length(pv.pc2.with_pos$calibrated.p))
range(pv.pc2.with_pos$FDR.p, na.rm = T)
pv.pc2.with_pos <- pv.pc2.with_pos %>% mutate(FDR.sig = ifelse(FDR.p < 0.05, "sig", "background"))
pv.pc2.with_pos %>% group_by(FDR.sig) %>% tally()

#3800 at 0.05

####################
#ggplots

#FDR with auto corrected P (algorithm GIF)
p1 = ggplot(pv.pc2.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = FDR.sig)) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "P value") +
theme(
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)
p1



png("figures/GxE/LFMM/bioclim.precip.PC1.png", width = 1080, height = 240)
p1
dev.off()

#######################
#PC3

#variable for test
X = sample_metadata.site_info$bioclim.PC3

#LFMM ridge
#mod.lfmm = lfmm_ridge(Y = Y, X = X, K = 2)
mod.lfmm.pc3 = lfmm_ridge(Y = Y, X = X, K = 2)

pv.pc3 <- lfmm_test(Y = Y,
X = X,
lfmm = mod.lfmm.pc3,
calibrate = "gif")


#Computing genomic inflation factor (GIF) based on calibrated z-scores (http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html) and Francois et al. 2016
lambda = median(pv.pc3$score^2, na.rm = T)/0.456
lambda #1.143442
adj.p.values = pchisq(pv.pc3$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
#Note that these calibrated scores are similar as pv$calibrated.pvalue
hist(pv.pc3$calibrated.pvalue)
hist(pv.pc3$pvalue)
#IN THIS CASE THE CALCULATED VALUES looks terrible and very conservative

#Try higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.pc3$score^2/1.35, df = 1, lower = FALSE)
hist(adj.p.values) #very conservative

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.pc3$score^2/1.05, df = 1, lower = FALSE)
hist(adj.p.values) #This looks good

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.pc3$score^2/0.85, df = 1, lower = FALSE)
hist(adj.p.values) #This looks good, but less conservative

#THIS IS SHOWING THAT THE GIF CALIBRATION IN THE ALGORITHM IS MORE CONSERVATIVE THAN LOWER VALUES OF LAMBDA
#However, the lower values have a correct distribution under null model
#Try at GIF = 0.85
adj.p.values = pchisq(pv.pc3$score^2/0.85, df = 1, lower = FALSE)
hist(adj.p.values)


#Join with actual positiion and chromosome
pv.pc3.with_pos = data.frame(calibrated.p = pv.pc3$calibrated.pvalue, effect_size = pv.pc3$B, SNP_pos, man.adj.p = adj.p.values) %>% left_join(., scf_lens)

#############################
#Outlier tests
# identify the 95% percentile
my_threshold <- quantile((pv.pc3.with_pos %>% filter(length > 100000))$calibrated.p, 0.025, na.rm = T)
# make an outlier column in the data.frame
pv.pc3.with_pos <- pv.pc3.with_pos %>% mutate(outlier = ifelse(calibrated.p < my_threshold, "outlier", "background"))
#Number of outliers
pv.pc3.with_pos %>% group_by(outlier) %>% tally()


#FDR correction
pv.pc3.with_pos$FDR.p = p.adjust(pv.pc3.with_pos$calibrated.p, method = "fdr", n = length(pv.pc3.with_pos$calibrated.p))
pv.pc3.with_pos <- pv.pc3.with_pos %>% mutate(FDR.sig = ifelse(FDR.p < 0.01, "sig", "background"))
pv.pc3.with_pos %>% group_by(FDR.sig) %>% tally()

#474

#FDR correction manual adjustment (lambda = 0.85)
pv.pc3.with_pos$FDR.p.man = p.adjust(pv.pc3.with_pos$man.adj.p, method = "fdr", n = length(pv.pc3.with_pos$man.adj.p))
pv.pc3.with_pos <- pv.pc3.with_pos %>% mutate(FDR.sig.man = ifelse(FDR.p.man < 0.05, "sig", "background"))
pv.pc3.with_pos %>% group_by(FDR.sig.man) %>% tally()

#506

####################
#ggplots


#FDR with auto corrected P (algorithm GIF)
p1 = ggplot(pv.ppt.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = FDR.sig)) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
#scale_x_continuous(labels = fancy_scientific, breaks = c(1, seq(from = 10^6, to = 6*10^6, by = 10^6)) ) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "P value") +
theme(
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)
p1

p2 = ggplot(pv.ppt.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man)) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
#scale_x_continuous(labels = fancy_scientific, breaks = c(1, seq(from = 10^6, to = 6*10^6, by = 10^6)) ) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = compose_trans("log10", "reverse"),
    labels = label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "P value") +
theme(
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)
p2

pdf("figures/LFMM.ppt.pdf", width = 18, height = 4)
p1
dev.off()


#######################
#Inection duration

#variable for test
X = sample_metadata.site_info$bioclim.PC4

#LFMM ridge
#mod.lfmm = lfmm_ridge(Y = Y, X = X, K = 2)
mod.lfmm.pc4 = lfmm_ridge(Y = Y, X = X, K = 2)

pv.pc4 <- lfmm_test(Y = Y,
X = X,
lfmm = mod.lfmm.pc4,
calibrate = "gif")

#Computing genomic inflation factor (GIF) based on calibrated z-scores (http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html) and Francois et al. 2016
lambda = median(pv.pc4$score^2, na.rm = T)/0.456
lambda #1.089663
adj.p.values = pchisq(pv.pc4$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
#Note that these calibrated scores are similar as pv$calibrated.pvalue
hist(pv.pc4$calibrated.pvalue)

#IN THIS CASE THE CALCULATED VALUES looks conservative

#Try higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.pc4$score^2/1.25, df = 1, lower = FALSE)
hist(adj.p.values) #very conservative

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.pc4$score^2/1, df = 0.9, lower = FALSE)
hist(adj.p.values) #This looks good

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.pc4$score^2/0.85, df = 1, lower = FALSE)
hist(adj.p.values)

#THIS IS SHOWING THAT THE GIF CALIBRATION IN THE ALGORITHM IS MORE CONSERVATIVE THAN LOWER VALUES OF LAMBDA
#However, the lower values have a correct distribution under null model
#Try at GIF = 0.1
adj.p.values = pchisq(pv.pc4$score^2/1, df = 1, lower = FALSE)
hist(adj.p.values)

#Read scaffold lengths
scf_lens = read.table("data/Nf_post_SPANDx/scaffold_lengths.csv", sep = ",", header = F)
colnames(scf_lens) = c("scaffold", "length")

#Join with actual positiion and chromosome
pv.pc4.with_pos = data.frame(calibrated.p = pv.pc4$calibrated.pvalue, effect_size = pv.pc4$B, SNP_pos, man.adj.p = adj.p.values) %>% left_join(., scf_lens)

#############################
#Outlier tests
# identify the 95% percentile
my_threshold <- quantile((pv.pc4.with_pos %>% filter(length > 100000))$calibrated.p, 0.025, na.rm = T)
# make an outlier column in the data.frame
pv.pc4.with_pos <- pv.pc4.with_pos %>% mutate(outlier = ifelse(calibrated.p < my_threshold, "outlier", "background"))
#Number of outliers
pv.pc4.with_pos %>% group_by(outlier) %>% tally()
#Example plot

#3223

#FDR correction
pv.pc4.with_pos$FDR.p = p.adjust(pv.pc4.with_pos$calibrated.p, method = "fdr", n = length(pv.pc4.with_pos$calibrated.p))
pv.pc4.with_pos <- pv.pc4.with_pos %>% mutate(FDR.sig = ifelse(FDR.p < 0.1, "sig", "background"))
pv.pc4.with_pos %>% group_by(FDR.sig) %>% tally()

#23

#FDR correction manual adjustment (lambda = 1)
pv.pc4.with_pos$FDR.p.man = p.adjust(pv.pc4.with_pos$man.adj.p, method = "fdr", n = length(pv.pc4.with_pos$man.adj.p))
pv.pc4.with_pos <- pv.pc4.with_pos %>% mutate(FDR.sig.man = ifelse(FDR.p.man < 0.05, "sig", "background"))
pv.pc4.with_pos %>% group_by(FDR.sig.man) %>% tally()

#23

####################
#ggplots

#Basic plot
ggplot(pv.pc4.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p)) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
#scale_x_continuous(labels = fancy_scientific, breaks = c(1, seq(from = 10^6, to = 6*10^6, by = 10^6)) ) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(trans = reverselog_trans(10), labels = trans_format('log10',math_format(10^.x))) +
my_gg_theme +
labs(x = "Position (Mbp)", y = "P value") +
theme(
strip.text.x = element_blank(),
axis.text.x = element_text(size = 8)
)

#Colored by outliers
ggplot(pv.pc4.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = outlier)) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
#scale_x_continuous(labels = fancy_scientific, breaks = c(1, seq(from = 10^6, to = 6*10^6, by = 10^6)) ) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(trans = reverselog_trans(10), labels = trans_format('log10',math_format(10^.x)), guide = "none") +
scale_color_manual(values = c("grey", "black")) +
my_gg_theme +
labs(x = "Position (Mbp)", y = "P value") +
theme(
strip.text.x = element_blank(),
axis.text.x = element_text(size = 8)
)

#FDR with auto corrected P (algorithm GIF)
p1 = ggplot(pv.pc4.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = FDR.sig)) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
#scale_x_continuous(labels = fancy_scientific, breaks = c(1, seq(from = 10^6, to = 6*10^6, by = 10^6)) ) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(trans = reverselog_trans(10), labels = trans_format('log10',math_format(10^.x))) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "P value") +
theme(
strip.text.x = element_blank(),
axis.text.x = element_text(size = 8)
#axis.text.x = element_blank(),
#axis.title.x = element_blank()
)

p2 = ggplot(pv.pc4.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man)) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
#scale_x_continuous(labels = fancy_scientific, breaks = c(1, seq(from = 10^6, to = 6*10^6, by = 10^6)) ) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(trans = reverselog_trans(10), labels = trans_format('log10',math_format(10^.x))) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "") +
theme(
strip.text.x = element_blank(),
#axis.text.x = element_text(size = 8)
axis.text.x = element_blank(),
axis.title.x = element_blank()
)

pdf("figures/LFMM.dur_inf.pdf", width = 18, height = 4)
p1
p2
dev.off()




#######################
#write tables

write.table(pv.hdd4.with_pos, "data/Nf_LFMM_tables/hdd4_lfmm.txt", quote = F, row.names = F, sep = "\t")
write.table(pv.ft.with_pos, "data/Nf_LFMM_tables/freezeThaw_lfmm.txt", quote = F, row.names = F, sep = "\t")
write.table(pv.ppt.with_pos, "data/Nf_LFMM_tables/ppt_lfmm.txt", quote = F, row.names = F, sep = "\t")
write.table(pv.dur_inf.with_pos, "data/Nf_LFMM_tables/dur_inf_lfmm.txt", quote = F, row.names = F, sep = "\t")



#Nice aligned plot of all three
library(gtable)
library(gridExtra)
library(grid)

############################################
##READ THE ABOVE TABLES BACK IN TO RUN BELOW

p1 = ggplot(pv.hdd4.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = FDR.sig)) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
#scale_x_continuous(labels = fancy_scientific, breaks = c(1, seq(from = 10^6, to = 6*10^6, by = 10^6)) ) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(trans = reverselog_trans(10), labels = trans_format('log10',math_format(10^.x))) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "P value") +
theme(
strip.text.x = element_blank(),
axis.text.x = element_text(size = 8)
#axis.text.x = element_blank(),
#axis.title.x = element_blank()
)
p1

p4 = ggplot(pv.hdd4.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man)) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
#scale_x_continuous(labels = fancy_scientific, breaks = c(1, seq(from = 10^6, to = 6*10^6, by = 10^6)) ) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(trans = reverselog_trans(10), labels = trans_format('log10',math_format(10^.x))) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "") +
theme(
strip.text.x = element_blank(),
#axis.text.x = element_text(size = 8)
axis.text.x = element_blank(),
axis.title.x = element_blank()
)

p2 = ggplot(pv.ft.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = FDR.sig)) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
#scale_x_continuous(labels = fancy_scientific, breaks = c(1, seq(from = 10^6, to = 6*10^6, by = 10^6)) ) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(trans = reverselog_trans(10), labels = trans_format('log10',math_format(10^.x))) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "P value") +
theme(
strip.text.x = element_blank(),
axis.text.x = element_text(size = 8)
#axis.text.x = element_blank(),
#axis.title.x = element_blank()
)

p5 = ggplot(pv.ft.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man)) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
#scale_x_continuous(labels = fancy_scientific, breaks = c(1, seq(from = 10^6, to = 6*10^6, by = 10^6)) ) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(trans = reverselog_trans(10), labels = trans_format('log10',math_format(10^.x))) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "") +
theme(
strip.text.x = element_blank(),
#axis.text.x = element_text(size = 8)
axis.text.x = element_blank(),
axis.title.x = element_blank()
)

p3 = ggplot(pv.ppt.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = FDR.sig)) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
#scale_x_continuous(labels = fancy_scientific, breaks = c(1, seq(from = 10^6, to = 6*10^6, by = 10^6)) ) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(trans = reverselog_trans(10), labels = trans_format('log10',math_format(10^.x))) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "P value") +
theme(
strip.text.x = element_blank(),
axis.text.x = element_text(size = 8)
#axis.text.x = element_blank(),
#axis.title.x = element_blank()
)

p6 = ggplot(pv.ppt.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man)) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
#scale_x_continuous(labels = fancy_scientific, breaks = c(1, seq(from = 10^6, to = 6*10^6, by = 10^6)) ) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(trans = reverselog_trans(10), labels = trans_format('log10',math_format(10^.x))) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "") +
theme(
strip.text.x = element_blank(),
#axis.text.x = element_text(size = 8)
axis.text.x = element_blank(),
axis.title.x = element_blank()
)

p7 = ggplot(pv.dur_inf.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = FDR.sig)) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
#scale_x_continuous(labels = fancy_scientific, breaks = c(1, seq(from = 10^6, to = 6*10^6, by = 10^6)) ) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(trans = reverselog_trans(10), labels = trans_format('log10',math_format(10^.x))) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "P value") +
theme(
strip.text.x = element_blank(),
axis.text.x = element_text(size = 8)
#axis.text.x = element_blank(),
#axis.title.x = element_blank()
)

p8 = ggplot(pv.dur_inf.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man)) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
#scale_x_continuous(labels = fancy_scientific, breaks = c(1, seq(from = 10^6, to = 6*10^6, by = 10^6)) ) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(trans = reverselog_trans(10), labels = trans_format('log10',math_format(10^.x))) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "") +
theme(
strip.text.x = element_blank(),
#axis.text.x = element_text(size = 8)
axis.text.x = element_blank(),
axis.title.x = element_blank()
)

########################
#AUTOMATIC ACLIPBRATED P-VALUES
plots = list(p1, p2, p3)

grobs <- lapply(plots, ggplotGrob)
# for gridExtra < v2.3, use do.call(gridExtra::rbind.gtable, grobs)
# for gridExtra >= v2.3 use:
g <- do.call(gridExtra::gtable_rbind, grobs)

#set relative panel heights
panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels] <- unit(c(1,1,1), "null")


grid.newpage()
grid.draw(g)


pdf("figures/LFMM.nongrowing_GDD.freeze_thaw.ppt.AUTO_GIF.pdf", width = 16, height = 6)
grid.newpage()
grid.draw(g)
dev.off()

#################################
#manual calibrated p-values
#################################

plots = list(p4, p5, p6)

grobs <- lapply(plots, ggplotGrob)
# for gridExtra < v2.3, use do.call(gridExtra::rbind.gtable, grobs)
# for gridExtra >= v2.3 use:
g <- do.call(gridExtra::gtable_rbind, grobs)

#set relative panel heights
panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels] <- unit(c(1,1,1), "null")


#grid.newpage()
#grid.draw(g)


pdf("figures/LFMM.nongrowing_GDD.freeze_thaw.ppt.MANUAL_CALIBRATION.pdf", width = 16, height = 6)
grid.newpage()
grid.draw(g)
dev.off()



##########################
##########################
##########################

#Find high density peaks in scf3

ggplot(pv.with_pos %>% filter(scaffold == "tig00000025_pilon" & position > 2.5*10^6 & position < 3.1*10^6), aes(x = position, y = calibrated.p, color = FDR.sig)) +
geom_point(alpha = 0.5, size = 1) +
scale_y_continuous(trans = reverselog_trans(10), labels = trans_format('log10',math_format(10^.x))) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "P value") +
theme(
strip.text.x = element_blank(),
axis.text.x = element_text(size = 8)
#axis.text.x = element_blank(),
#axis.title.x = element_blank()
)



#Should find a way to do this programmatically, but...

ggplot(pv.with_pos %>% filter(scaffold == "tig00000025_pilon" & position > 2.5*10^6 & position < 2.6*10^6), aes(x = position, y = calibrated.p, color = FDR.sig)) +
geom_point(alpha = 0.5, size = 1) +
scale_y_continuous(trans = reverselog_trans(10), labels = trans_format('log10',math_format(10^.x))) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "P value") +
theme(
strip.text.x = element_blank(),
axis.text.x = element_text(size = 8)
#axis.text.x = element_blank(),
#axis.title.x = element_blank()
)

pv.with_pos %>% filter(scaffold == "tig00000025_pilon" & position > 2.5*10^6 & position < 2.6*10^6 & FDR.sig == "sig") %>% select(position) %>% nrow
#39
pv.with_pos %>% filter(scaffold == "tig00000025_pilon" & position > 2.5*10^6 & position < 2.6*10^6 & FDR.sig == "sig") %>% select(position) %>% range
#2524273-2544922


ggplot(pv.with_pos %>% filter(scaffold == "tig00000025_pilon" & position > 3.05*10^6 & position < 3.1*10^6), aes(x = position, y = calibrated.p, color = FDR.sig)) +
geom_point(alpha = 0.5, size = 1) +
scale_y_continuous(trans = reverselog_trans(10), labels = trans_format('log10',math_format(10^.x))) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "P value") +
theme(
strip.text.x = element_blank(),
axis.text.x = element_text(size = 8)
#axis.text.x = element_blank(),
#axis.title.x = element_blank()
)

pv.with_pos %>% filter(scaffold == "tig00000025_pilon" & position > 3.05*10^6 & position < 3.1*10^6 & FDR.sig == "sig") %>% select(position) %>% nrow
#98
pv.with_pos %>% filter(scaffold == "tig00000025_pilon" & position > 3.05*10^6 & position < 3.1*10^6 & FDR.sig == "sig") %>% select(position) %>% range
#3070669-3078240
