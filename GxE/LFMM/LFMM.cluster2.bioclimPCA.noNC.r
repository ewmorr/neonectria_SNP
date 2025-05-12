library(lfmm)
library(dplyr)
library(ggplot2)
library(scales)
source("library/ggplot_theme.txt")

#read site metadata
site.info = read.csv("data/sample_metadata/site_info.csv")
site.bioclim = read.csv("data/sample_metadata/Nf.site_bioclim.csv")
site.bioclim = site.bioclim %>% filter(state != "WV" & state != "NC")
site.bioclim.PC = read.csv("data/sample_metadata/Nf.sites_bioclim_PC.no_WV_no_NC.csv", header = T)
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
sample_ids = read.table("data/Nf/final_tables/rm_dups/FINAL_snp.gwas_analyses.cluster2_no_NC.sampleIDs", header = F)
sample_metadata.Nf.sorted = sample_metadata.Nf[sample_ids$V1,]
nrow(sample_metadata.Nf.sorted)
sample_metadata.site_info = left_join(sample_metadata.Nf.sorted, site_metadata)
is.na(sample_metadata.site_info$bio1) %>% sum()

row_ids = which(sample_ids$V1 %in% sample_metadata.site_info$Sequence_label)

#The genotype data can simply be read in as a matrix (according the docs)
#OR can try loading LEA and using readLfmm()
Y = as.matrix(read.table("data/Nf/final_tables/rm_dups/FINAL_snp.gwas_analyses.cluster2_no_NC.lfmm", header = F))
Y.filtered = Y[row_ids,]
nrow(Y.filtered)
ncol(Y.filtered)

SNP_pos = read.table("data/Nf/final_tables/rm_dups/FINAL_snp.gwas_analyses.cluster2_no_NC.recode.map")
nrow(SNP_pos)
SNP_pos = SNP_pos[c(1,4)]
colnames(SNP_pos) = c("scaffold", "position")


###################
###################
###################
#MAC filter
ref_sum = colSums(Y.filtered == 0)
alt_sum = colSums(Y.filtered == 1)

#filter for MAC ge 3
nrow(Y.filtered)
minMAC = 3

sum(ref_sum < minMAC)
# 112892
sum(alt_sum < minMAC)
ncol(Y.filtered)
# 424811
ncol(Y.filtered) - sum(ref_sum < minMAC)
# 311919

which(ref_sum < minMAC) 
which(alt_sum < minMAC)# there are none
rm_cols = which(ref_sum < minMAC)

length(rm_cols)
#112892
ncol(Y.filtered)
length(rm_cols)/ncol(Y.filtered)
# 0.2657464
ncol(Y.filtered)-length(rm_cols)
# 311919
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

nrow(sample_metadata.site_info)
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
lambda #1.077236
adj.p.values = pchisq(pv$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
#Note that these calibrated scores are similar as pv$calibrated.pvalue
hist(pv$calibrated.pvalue) #this looks basically perfect
hist(pv$pvalue)
#Try higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv$score^2/1.15, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv$score^2/0.95, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv$score^2/1.0, df = 1, lower = FALSE)
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
#7798
#FDR correction
#This is based on the auto calibartion
pv.with_pos$FDR.p = p.adjust(pv.with_pos$calibrated.p, method = "fdr", n = length(pv.with_pos$calibrated.p))
pv.with_pos <- pv.with_pos %>% mutate(FDR.sig = ifelse(FDR.p < 0.05, "sig", "background"))
pv.with_pos %>% group_by(FDR.sig) %>% tally()

#0  SNPs identified as significant after FDR correction


#FDR correction
#This is based on the manual GIF adjustment
pv.with_pos$FDR.p.man = p.adjust(pv.with_pos$man.adj.p, method = "fdr", n = length(pv.with_pos$man.adj.p))
pv.with_pos <- pv.with_pos %>% mutate(FDR.sig.man = ifelse(FDR.p.man < 0.05, "sig", "background"))
pv.with_pos %>% group_by(FDR.sig.man) %>% tally()


pv.PC1.with_pos = pv.with_pos
#0 SNPs identified as significant after FDR correction


pv.PC1.with_pos %>% head()
pv.PC1.with_pos %>% filter(outlier == "outlier") %>% nrow()
#7798 
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

png("figures/GxE/LFMM/no_NC.bioclim.temp.PC1.png", width = 1080, height = 240)
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
lambda #1.195644
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
adj.p.values = pchisq(pv.pc2$score^2/1.05, df = 1, lower = FALSE)
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

#8 at 0.05

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



png("figures/GxE/LFMM/no_NC.bioclim.precip.PC1.png", width = 1080, height = 240)
p1
dev.off()

#######################
#Precip pc2

#variable for test
X = sample_metadata.site_info$precip.PC2

#LFMM ridge
#mod.lfmm = lfmm_ridge(Y = Y, X = X, K = 2)
mod.lfmm.precip.PC2 = lfmm_ridge(Y = Y.filtered, X = X, K = 2)

pv.precip.PC2 <- lfmm_test(Y = Y.filtered,
X = X,
lfmm = mod.lfmm.precip.PC2,
calibrate = "gif")


#Computing genomic inflation factor (GIF) based on calibrated z-scores (http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html) and Francois et al. 2016
lambda = median(pv.precip.PC2$score^2, na.rm = T)/0.456
lambda #0.9586345
adj.p.values = pchisq(pv.precip.PC2$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
#Note that these calibrated scores are similar as pv$calibrated.pvalue
hist(pv.precip.PC2$calibrated.pvalue)
hist(pv.precip.PC2$pvalue)
#IN THIS CASE THE CALCULATED VALUES looks ok

#Try higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.precip.PC2$score^2/1.15, df = 1, lower = FALSE)
hist(adj.p.values) #very conservative

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.precip.PC2$score^2/1.05, df = 1, lower = FALSE)
hist(adj.p.values) #This looks conservative

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.precip.PC2$score^2/0.85, df = 1, lower = FALSE)
hist(adj.p.values) #This looks good, but less conservative

#THIS IS SHOWING THAT THE GIF CALIBRATION IN THE ALGORITHM IS MORE CONSERVATIVE THAN LOWER VALUES OF LAMBDA
#However, the lower values have a correct distribution under null model
#Try at GIF = 1.0
adj.p.values = pchisq(pv.precip.PC2$score^2/0.90, df = 1, lower = FALSE)
hist(adj.p.values)


#Join with actual positiion and chromosome
nrow(SNP_pos)
nrow(pv.precip.PC2$calibrated.pvalue)
pv.precip.PC2.with_pos = data.frame(calibrated.p = pv.precip.PC2$calibrated.pvalue, effect_size = pv.precip.PC2$B, SNP_pos, man.adj.p = adj.p.values) %>% left_join(., scf_lens)

#############################
#Outlier tests
# identify the 95% percentile
my_threshold <- quantile((pv.precip.PC2.with_pos %>% filter(length > 100000))$calibrated.p, 0.025, na.rm = T)
# make an outlier column in the data.frame
pv.precip.PC2.with_pos <- pv.precip.PC2.with_pos %>% mutate(outlier = ifelse(calibrated.p < my_threshold, "outlier", "background"))
#Number of outliers
pv.precip.PC2.with_pos %>% group_by(outlier) %>% tally()


#FDR correction
pv.precip.PC2.with_pos$FDR.p = p.adjust(pv.precip.PC2.with_pos$calibrated.p, method = "fdr", n = length(pv.precip.PC2.with_pos$calibrated.p))
pv.precip.PC2.with_pos <- pv.precip.PC2.with_pos %>% mutate(FDR.sig = ifelse(FDR.p < 0.05, "sig", "background"))
pv.precip.PC2.with_pos %>% group_by(FDR.sig) %>% tally()

#none

#FDR correction manual adjustment (lambda = 1)
pv.precip.PC2.with_pos$FDR.p.man = p.adjust(pv.precip.PC2.with_pos$man.adj.p, method = "fdr", n = length(pv.precip.PC2.with_pos$man.adj.p))
pv.precip.PC2.with_pos <- pv.precip.PC2.with_pos %>% mutate(FDR.sig.man = ifelse(FDR.p.man < 0.05, "sig", "background"))
pv.precip.PC2.with_pos %>% group_by(FDR.sig.man) %>% tally()

#1

####################
#ggplots


#FDR with auto corrected P (algorithm GIF)
p1 = ggplot(pv.precip.PC2.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = FDR.sig)) +
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

p2 = ggplot(pv.precip.PC2.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man)) +
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

png("figures/GxE/LFMM/no_NC.bioclim.precip.PC2.png", width = 1080, height = 240)
p1
dev.off()


#######################
#bio2

#variable for test
X = sample_metadata.site_info$bio2

#LFMM ridge
#mod.lfmm = lfmm_ridge(Y = Y, X = X, K = 2)
mod.lfmm.bio2 = lfmm_ridge(Y = Y.filtered, X = X, K = 2)

pv.bio2 <- lfmm_test(Y = Y.filtered,
X = X,
lfmm = mod.lfmm.bio2,
calibrate = "gif")

#Computing genomic inflation factor (GIF) based on calibrated z-scores (http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html) and Francois et al. 2016
lambda = median(pv.bio2$score^2, na.rm = T)/0.456
lambda #1.043716
adj.p.values = pchisq(pv.bio2$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
#Note that these calibrated scores are similar as pv$calibrated.pvalue
hist(pv.bio2$calibrated.pvalue)

#IN THIS CASE THE CALCULATED VALUES looks conservative

#Try higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.bio2$score^2/1.25, df = 1, lower = FALSE)
hist(adj.p.values) #very conservative

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.bio2$score^2/1.15, df = 1, lower = FALSE)
hist(adj.p.values) #This looks very conservative


#THIS IS SHOWING THAT THE GIF CALIBRATION IN THE ALGORITHM IS MORE CONSERVATIVE THAN LOWER VALUES OF LAMBDA
#However, the lower values have a correct distribution under null model
#Try at GIF = 1.25
adj.p.values = pchisq(pv.bio2$score^2/0.95, df = 1, lower = FALSE)
hist(adj.p.values)

#Join with actual positiion and chromosome
pv.bio2.with_pos = data.frame(calibrated.p = pv.bio2$calibrated.pvalue, effect_size = pv.bio2$B, SNP_pos, man.adj.p = adj.p.values) %>% left_join(., scf_lens)

#############################
#Outlier tests
# identify the 95% percentile
my_threshold <- quantile((pv.bio2.with_pos %>% filter(length > 100000))$calibrated.p, 0.025, na.rm = T)
# make an outlier column in the data.frame
pv.bio2.with_pos <- pv.bio2.with_pos %>% mutate(outlier = ifelse(calibrated.p < my_threshold, "outlier", "background"))
#Number of outliers
pv.bio2.with_pos %>% group_by(outlier) %>% tally()
#Example plot

#FDR correction
pv.bio2.with_pos$FDR.p = p.adjust(pv.bio2.with_pos$calibrated.p, method = "fdr", n = length(pv.bio2.with_pos$calibrated.p))
pv.bio2.with_pos <- pv.bio2.with_pos %>% mutate(FDR.sig = ifelse(FDR.p < 0.05, "sig", "background"))
pv.bio2.with_pos %>% group_by(FDR.sig) %>% tally()

#none

#FDR correction manual adjustment (lambda = 1)
pv.bio2.with_pos$FDR.p.man = p.adjust(pv.bio2.with_pos$man.adj.p, method = "fdr", n = length(pv.bio2.with_pos$man.adj.p))
pv.bio2.with_pos <- pv.bio2.with_pos %>% mutate(FDR.sig.man = ifelse(FDR.p.man < 0.05, "sig", "background"))
pv.bio2.with_pos %>% group_by(FDR.sig.man) %>% tally()

#0

####################
#ggplots

#FDR with auto corrected P (algorithm GIF)
p1 = ggplot(pv.bio2.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = FDR.sig)) +
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

p2 = ggplot(pv.bio2.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man)) +
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
labs(x = "Position (Mbp)", y = "") +
theme(
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)

png("figures/GxE/LFMM/no_NCbioclim.bio2.png", width = 1080, height = 240)
p1
dev.off()


#######################
#duration_infection

#variable for test
X = sample_metadata.site_info$duration_infection

#LFMM ridge
#mod.lfmm = lfmm_ridge(Y = Y, X = X, K = 2)
mod.lfmm.duration_infection = lfmm_ridge(Y = Y.filtered, X = X, K = 2)

pv.duration_infection <- lfmm_test(Y = Y.filtered,
X = X,
lfmm = mod.lfmm.duration_infection,
calibrate = "gif")

#Computing genomic inflation factor (GIF) based on calibrated z-scores (http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html) and Francois et al. 2016
lambda = median(pv.duration_infection$score^2, na.rm = T)/0.456
lambda #0.990194
adj.p.values = pchisq(pv.duration_infection$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
#Note that these calibrated scores are similar as pv$calibrated.pvalue
hist(pv.duration_infection$calibrated.pvalue)

#IN THIS CASE THE CALCULATED VALUES looks conservative

#Try higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.duration_infection$score^2/1.05, df = 1, lower = FALSE)
hist(adj.p.values) #very conservative

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.duration_infection$score^2/0.85, df = 1, lower = FALSE)
hist(adj.p.values) #This looks good


#THIS IS SHOWING THAT THE GIF CALIBRATION IN THE ALGORITHM IS MORE CONSERVATIVE THAN LOWER VALUES OF LAMBDA
#However, the lower values have a correct distribution under null model
#Try at GIF = 1.05
adj.p.values = pchisq(pv.duration_infection$score^2/1.05, df = 1, lower = FALSE)
hist(adj.p.values)

#Join with actual positiion and chromosome
pv.duration_infection.with_pos = data.frame(calibrated.p = pv.duration_infection$calibrated.pvalue, effect_size = pv.duration_infection$B, SNP_pos, man.adj.p = adj.p.values) %>% left_join(., scf_lens)

#############################
#Outlier tests
# identify the 95% percentile
my_threshold <- quantile((pv.duration_infection.with_pos %>% filter(length > 100000))$calibrated.p, 0.025, na.rm = T)
# make an outlier column in the data.frame
pv.duration_infection.with_pos <- pv.duration_infection.with_pos %>% mutate(outlier = ifelse(calibrated.p < my_threshold, "outlier", "background"))
#Number of outliers
pv.duration_infection.with_pos %>% group_by(outlier) %>% tally()
#Example plot

#FDR correction
pv.duration_infection.with_pos$FDR.p = p.adjust(pv.duration_infection.with_pos$calibrated.p, method = "fdr", n = length(pv.duration_infection.with_pos$calibrated.p))
pv.duration_infection.with_pos <- pv.duration_infection.with_pos %>% mutate(FDR.sig = ifelse(FDR.p < 0.05, "sig", "background"))
pv.duration_infection.with_pos %>% group_by(FDR.sig) %>% tally()

#1

#FDR correction manual adjustment (lambda = 1)
pv.duration_infection.with_pos$FDR.p.man = p.adjust(pv.duration_infection.with_pos$man.adj.p, method = "fdr", n = length(pv.duration_infection.with_pos$man.adj.p))
pv.duration_infection.with_pos <- pv.duration_infection.with_pos %>% mutate(FDR.sig.man = ifelse(FDR.p.man < 0.05, "sig", "background"))
pv.duration_infection.with_pos %>% group_by(FDR.sig.man) %>% tally()

#6

####################
#ggplots

#FDR with auto corrected P (algorithm GIF)
p1 = ggplot(pv.duration_infection.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = FDR.sig)) +
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
    axis.text.x = element_text(size = 8)
)
p1

p2 = ggplot(pv.duration_infection.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man)) +
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
labs(x = "Position (Mbp)", y = "") +
theme(
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)
p2

png("figures/GxE/LFMM/no_NC.duration_infection.png", width = 1080, height = 240)
p1
dev.off()



#######################
#write tables

write.table(pv.PC1.with_pos, "data/Nf/GxE/LFMM/no_NC.temp.pc1.lfmm.txt", quote = F, row.names = F, sep = "\t")
write.table(pv.pc2.with_pos, "data/Nf/GxE/LFMM/no_NC.precip.pc1.lfmm.txt", quote = F, row.names = F, sep = "\t")
write.table(pv.precip.PC2.with_pos, "data/Nf/GxE/LFMM/no_NC.precip.pc2.lfmm.txt", quote = F, row.names = F, sep = "\t")
write.table(pv.bio2.with_pos, "data/Nf/GxE/LFMM/no_NC.bio2.lfmm.txt", quote = F, row.names = F, sep = "\t")
write.table(pv.duration_infection.with_pos, "data/Nf/GxE/LFMM/no_NC.dur_inf.lfmm.txt", quote = F, row.names = F, sep = "\t")



#Nice aligned plot of all three
library(gtable)
library(gridExtra)
library(grid)

############################################
##READ THE ABOVE TABLES BACK IN TO RUN BELOW

p1 = ggplot(pv.PC1.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = FDR.sig)) +
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

p2 = ggplot(pv.pc2.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = FDR.sig)) +
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

p3 = ggplot(pv.precip.PC2.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = FDR.sig)) +
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

p4 = ggplot(pv.bio2.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = FDR.sig)) +
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

p5 = ggplot(pv.duration_infection.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = FDR.sig)) +
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
    axis.text.x = element_text(size = 8)
)

########################
#AUTOMATIC ACLIPBRATED P-VALUES
plots = list(p1, p2, p3, p4, p5)

grobs <- lapply(plots, ggplotGrob)
# for gridExtra < v2.3, use do.call(gridExtra::rbind.gtable, grobs)
# for gridExtra >= v2.3 use:
g <- do.call(gridExtra::gtable_rbind, grobs)

#set relative panel heights
panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels] <- unit(rep(1,5), "null")


grid.newpage()
grid.draw(g)


png("figures/GxE/LFMM/no_NC.cluster2.tempPC1-precipPC1-precipPC2-bio2-durInf.AUTO_GIF.png", width = 1080, height = 640)
grid.newpage()
grid.draw(g)
dev.off()

