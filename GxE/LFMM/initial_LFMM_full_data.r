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
site.info = site.info %>% filter(country == "USA")
site.GDD = read.table("data/PRISM_dailies/site_climate.GDD.txt", header = T)
site.climate = read.csv("data/sample_metadata/sites_climate.csv", header = T)
site.GDD$freezeThaw.annual_mean = site.GDD$freezeThaw.mean_growing + site.GDD$freezeThaw.mean_nongrowing

site_metadata = left_join(site.GDD, site.info %>% select(Site, lat, lon, duration_infection), by = "Site") %>%
    left_join(., site.climate %>% select(Site, tmin, tmax, ppt, MAT, lat, lon, state), by = c("Site", "lat", "lon") ) %>%
    select(-Site) %>%
    unique()
site_metadata = site_metadata[-16,] #remove extra nh.ccm
#sample metadata. Needs to be sorted by PED .sampleIDs info
sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
rownames(sample_metadata.Nf) = sample_metadata.Nf$Sequence_label
#PED sample IDs
sample_ids = read.table("data/Nf/final_tables/rm_dups/FINAL_snp.mac_ge2.biallele.gwas_analyses.sampleIDs", header = F)
sample_metadata.Nf.sorted = sample_metadata.Nf[sample_ids$V1,]
sample_metadata.Nf.US = sample_metadata.Nf.sorted %>% filter(country == "USA")
sample_metadata.site_info = left_join(sample_metadata.Nf.US, site_metadata)
is.na(sample_metadata.site_info$HDD4.mean_nongrowing) %>% sum()

row_ids = which(sample_ids$V1 %in% sample_metadata.Nf.US$Sequence_label)

#The genotype data can simply be read in as a matrix (according the docs)
#OR can try loading LEA and using readLfmm()
Y = as.matrix(read.table("data/Nf/final_tables/rm_dups/FINAL_snp.mac_ge2.biallele.gwas_analyses.lfmm", header = F))
Y.filtered = Y[row_ids,]

SNP_pos = read.table("data/Nf/final_tables/rm_dups/FINAL_snp.mac_ge2.biallele.gwas_analyses.recode.map")
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
minMAC = 10

sum(ref_sum < minMAC)
# 161795
sum(alt_sum < minMAC)
ncol(Y.filtered)
# 424811
ncol(Y.filtered) - sum(ref_sum < minMAC)
# 263016

which(ref_sum < minMAC) 
which(alt_sum < minMAC)# there are none
rm_cols = which(ref_sum < minMAC)

length(rm_cols)
ncol(Y.filtered)
length(rm_cols)/ncol(Y.filtered)
ncol(Y.filtered)-length(rm_cols)
# 263016
Y.filteredMAC = Y.filtered[,-rm_cols]
#filter for greater than 25% missing data
sum(colSums(Y.filteredMAC == 9)/nrow(Y.filteredMAC) > 0.05)
rm_cols = which(colSums(Y.filteredMAC == 9)/nrow(Y.filteredMAC) > 0.05)
Y.filteredNA = Y.filteredMAC[,-rm_cols]

ncol(Y.filteredNA)
#207548
str(Y.filteredNA)

plot(colSums(Y == 9)/nrow(Y))
plot(colSums(Y == 0))
plot(colSums(Y == 1))
###################
###################
###################

#
#principal components analysis for K
pc = prcomp(Y.filtered)
plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
points(4,pc$sdev[7]^2, type = "h", lwd = 3, col = "blue")


#Y = Y.filteredNA
#######################
#NONGROWING SEASON HDD

#variable for test
X = (sample_metadata.site_info$HDD4.mean_nongrowing)

class(X)
class(Y.filtered)
sd(X)
mean(X)
#LFMM ridge

mod.lfmm = lfmm_ridge(Y = Y.filtered, X = X, K = 4, lambda = 1) #using K = 4 based on PCA (STRUCTURE suggests 2 but this method is PCA based)
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
xlab = "Probe", ylab = "-Log P",
col = "grey")

plot((pv$B),
pch = 19,
cex = .3,
xlab = "Probe", ylab = "Effect size",
col = "grey")

#Computing genomic inflation factor (GIF) based on calibrated z-scores (http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html) and Francois et al. 2016
lambda = median(pv$score^2)/0.456
lambda #1.509562
adj.p.values = pchisq(pv$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
#Note that these calibrated scores are similar as pv$calibrated.pvalue
hist(pv$calibrated.pvalue) #this is quite conservative
hist(pv$pvalue)
#Try higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv$score^2/1.25, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv$score^2/0.95, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv$score^2/0.85, df = 1, lower = FALSE)
hist(adj.p.values)

#THIS IS SHOWING THAT THE GIF CALIBRATION IN THE ALGORITHM IS MORE CONSERVATIVE THAN LOWER VALUES OF LAMBDA
#However, the lower values have a correct distribution under null model
#Try at GIF = 1.15
adj.p.values = pchisq(pv$score^2/1.15, df = 1, lower = FALSE)
hist(adj.p.values)

#Read scaffold lengths
scf_lens = read.table("data/Nf_SPANDx_all_seqs/scaffold_lengths.csv", sep = ",", header = F)
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
pv.with_pos %>% group_by(outlier) %>% tally()
#3216
#FDR correction
#This is based on the auto calibartion
pv.with_pos$FDR.p = p.adjust(pv.with_pos$calibrated.p, method = "fdr", n = length(pv.with_pos$calibrated.p))
pv.with_pos <- pv.with_pos %>% mutate(FDR.sig = ifelse(FDR.p < 0.1, "sig", "background"))
pv.with_pos %>% group_by(FDR.sig) %>% tally()

#178 of 128629 SNPs identified as significant after FDR correction

#FDR correction
#This is based on the manual GIF adjustment
pv.with_pos$FDR.p.man = p.adjust(pv.with_pos$man.adj.p, method = "fdr", n = length(pv.with_pos$man.adj.p))
pv.with_pos <- pv.with_pos %>% mutate(FDR.sig.man = ifelse(FDR.p.man < 0.05, "sig", "background"))
pv.with_pos %>% group_by(FDR.sig.man) %>% tally()
pv.hdd4.with_pos = pv.with_pos
#173 of 128629 SNPs identified as significant after FDR correction

####################
#ggplots

#Basic plot
ggplot(pv.hdd4.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p)) +
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
    #axis.text.x = element_text(angle = 85, size = 10, hjust = 1)
)

#Colored by outliers
ggplot(pv.hdd4.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = outlier)) +
#facet_wrap(~scaffold) +
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
#axis.text.x = element_text(angle = 85, size = 10, hjust = 1)
)

#FDR with auto corrected P (algorithm GIF)
p1 = ggplot(pv.hdd4.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = FDR.sig)) +
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

p2 = ggplot(pv.hdd4.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man)) +
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


#FDR correction maunally adjusted P (GIF = 0.95)
ggplot(pv.hdd4.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man)) +
#facet_wrap(~scaffold) +
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
#axis.text.x = element_text(angle = 85, size = 10, hjust = 1)
)

pdf("figures/LFMM.nongrowing_season_GDD.05312023.pdf", width = 18, height = 4)
p1
p2
dev.off()

#######################
#NONGROWING SEASON freeze-thaw
#######
#NOW SWITCHED TO MULTI ANNUAL MEAN OF TOTAL FREEZE THAW EVENTS

#variable for test
X = sample_metadata.site_info$freezeThaw.annual_mean

#LFMM ridge
#mod.lfmm = lfmm_ridge(Y = Y, X = X, K = 2)
mod.lfmm.ft = lfmm_ridge(Y = Y, X = X, K = 4)

pv.ft <- lfmm_test(Y = Y,
X = X,
lfmm = mod.lfmm.ft,
calibrate = "gif")

#Example plots
plot(-log10(pv.ft$calibrated.pvalue),
pch = 19,
cex = .3,
xlab = "Probe", ylab = "-Log P",
col = "grey")

plot((pv.ft$B),
pch = 19,
cex = .3,
xlab = "Probe", ylab = "Effect size",
col = "grey")

#Computing genomic inflation factor (GIF) based on calibrated z-scores (http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html) and Francois et al. 2016
lambda = median(pv.ft$score^2)/0.456
lambda #1.36
adj.p.values = pchisq(pv.ft$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
#Note that these calibrated scores are similar as pv$calibrated.pvalue
hist(pv.ft$calibrated.pvalue)
hist(pv.ft$pvalue)
#IN THIS CASE THE CALCULATED VALUES AREADY LOOK GOOD 

#Try higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.ft$score^2/1.45, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.ft$score^2/1.25, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.ft$score^2/0.15, df = 1, lower = FALSE)
hist(adj.p.values)

#THIS IS SHOWING THAT THE GIF CALIBRATION IN THE ALGORITHM IS Good
#However, the lower values have a correct distribution under null model
#Try at GIF = 1.25
adj.p.values = pchisq(pv.ft$score^2/1.35, df = 1, lower = FALSE)
hist(adj.p.values)

#Read scaffold lengths
scf_lens = read.table("data/Nf_post_SPANDx/scaffold_lengths.csv", sep = ",", header = F)
colnames(scf_lens) = c("scaffold", "length")

#Join with actual positiion and chromosome
pv.ft.with_pos = data.frame(calibrated.p = pv.ft$calibrated.pvalue, effect_size = pv.ft$B, SNP_pos, man.adj.p = adj.p.values) %>% left_join(., scf_lens)

#############################
#Outlier tests
# identify the 95% percentile
my_threshold <- quantile((pv.ft.with_pos)$calibrated.p, 0.025, na.rm = T)
# make an outlier column in the data.frame
pv.ft.with_pos <- pv.ft.with_pos %>% mutate(outlier = ifelse(calibrated.p < my_threshold, "outlier", "background"))
#Number of outliers
pv.ft.with_pos %>% group_by(outlier) %>% tally()

#3216

#FDR correction
pv.ft.with_pos$FDR.p = p.adjust(pv.ft.with_pos$calibrated.p, method = "fdr", n = length(pv.ft.with_pos$calibrated.p))
pv.ft.with_pos <- pv.ft.with_pos %>% mutate(FDR.sig = ifelse(FDR.p < 0.1, "sig", "background"))
pv.ft.with_pos %>% group_by(FDR.sig) %>% tally()

#65

#FDR correction manual adjustment (lambda = 1.15)
pv.ft.with_pos$FDR.p.man = p.adjust(pv.ft.with_pos$man.adj.p, method = "fdr", n = length(pv.ft.with_pos$man.adj.p))
pv.ft.with_pos <- pv.ft.with_pos %>% mutate(FDR.sig.man = ifelse(FDR.p.man < 0.05, "sig", "background"))
pv.ft.with_pos %>% group_by(FDR.sig.man) %>% tally()

#50

####################
#ggplots

#Basic plot
ggplot(pv.ft.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p)) +
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
#axis.text.x = element_text(angle = 85, size = 10, hjust = 1)
)

#Colored by outliers
ggplot(pv.ft.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = outlier)) +
#facet_wrap(~scaffold) +
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
#axis.text.x = element_text(angle = 85, size = 10, hjust = 1)
)

#FDR with auto corrected P (algorithm GIF)
p1 = ggplot(pv.ft.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = FDR.sig)) +
#facet_wrap(~scaffold) +
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
#axis.text.x = element_text(size = 8)
axis.text.x = element_blank(),
axis.title.x = element_blank()
)

p2 = ggplot(pv.ft.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man)) +
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


#FDR correction maunally adjusted P (GIF = 1.35)
ggplot(pv.ft.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man)) +
#facet_wrap(~scaffold) +
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
#axis.text.x = element_text(angle = 85, size = 10, hjust = 1)
)

pdf("figures/LFMM.nongrowing_season_freeze_thaw.01112023.pdf", width = 18, height = 4)
p1
p2
dev.off()

#######################
#Annual precip (ppt_

#variable for test
X = sample_metadata.site_info$ppt

#LFMM ridge
#mod.lfmm = lfmm_ridge(Y = Y, X = X, K = 2)
mod.lfmm.ppt = lfmm_ridge(Y = Y, X = X, K = 4)

pv.ppt <- lfmm_test(Y = Y,
X = X,
lfmm = mod.lfmm.ppt,
calibrate = "gif")

#Example plots
plot(-log10(pv.ppt$calibrated.pvalue),
pch = 19,
cex = .3,
xlab = "Probe", ylab = "-Log P",
col = "grey")

plot((pv.ppt$B),
pch = 19,
cex = .3,
xlab = "Probe", ylab = "Effect size",
col = "grey")

#Computing genomic inflation factor (GIF) based on calibrated z-scores (http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html) and Francois et al. 2016
lambda = median(pv.ppt$score^2)/0.456
lambda #0.95
adj.p.values = pchisq(pv.ppt$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
#Note that these calibrated scores are similar as pv$calibrated.pvalue
hist(pv.ppt$calibrated.pvalue)
hist(pv.ppt$pvalue)
#IN THIS CASE THE CALCULATED VALUES looks pretty good

#Try higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.ppt$score^2/1.15, df = 1, lower = FALSE)
hist(adj.p.values) #very conservative

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.ppt$score^2/0.95, df = 1, lower = FALSE)
hist(adj.p.values) #This looks good

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.ppt$score^2/0.85, df = 1, lower = FALSE)
hist(adj.p.values) #This looks good, but less conservative

#THIS IS SHOWING THAT THE GIF CALIBRATION IN THE ALGORITHM IS MORE CONSERVATIVE THAN LOWER VALUES OF LAMBDA
#However, the lower values have a correct distribution under null model
#Try at GIF = 0.85
adj.p.values = pchisq(pv.ppt$score^2/0.85, df = 1, lower = FALSE)
hist(adj.p.values)

#Read scaffold lengths
scf_lens = read.table("data/Nf_post_SPANDx/scaffold_lengths.csv", sep = ",", header = F)
colnames(scf_lens) = c("scaffold", "length")

#Join with actual positiion and chromosome
pv.ppt.with_pos = data.frame(calibrated.p = pv.ppt$calibrated.pvalue, effect_size = pv.ppt$B, SNP_pos, man.adj.p = adj.p.values) %>% left_join(., scf_lens)

#############################
#Outlier tests
# identify the 95% percentile
my_threshold <- quantile((pv.ppt.with_pos %>% filter(length > 100000))$calibrated.p, 0.025, na.rm = T)
# make an outlier column in the data.frame
pv.ppt.with_pos <- pv.ppt.with_pos %>% mutate(outlier = ifelse(calibrated.p < my_threshold, "outlier", "background"))
#Number of outliers
pv.ppt.with_pos %>% group_by(outlier) %>% tally()

#3214

#FDR correction
pv.ppt.with_pos$FDR.p = p.adjust(pv.ppt.with_pos$calibrated.p, method = "fdr", n = length(pv.ppt.with_pos$calibrated.p))
pv.ppt.with_pos <- pv.ppt.with_pos %>% mutate(FDR.sig = ifelse(FDR.p < 0.1, "sig", "background"))
pv.ppt.with_pos %>% group_by(FDR.sig) %>% tally()

#474

#FDR correction manual adjustment (lambda = 0.85)
pv.ppt.with_pos$FDR.p.man = p.adjust(pv.ppt.with_pos$man.adj.p, method = "fdr", n = length(pv.ppt.with_pos$man.adj.p))
pv.ppt.with_pos <- pv.ppt.with_pos %>% mutate(FDR.sig.man = ifelse(FDR.p.man < 0.05, "sig", "background"))
pv.ppt.with_pos %>% group_by(FDR.sig.man) %>% tally()

#506

####################
#ggplots

#Basic plot
ggplot(pv.ppt.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p)) +
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
#axis.text.x = element_text(angle = 85, size = 10, hjust = 1)
)

#Colored by outliers
ggplot(pv.ppt.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = outlier)) +
#facet_wrap(~scaffold) +
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
#axis.text.x = element_text(angle = 85, size = 10, hjust = 1)
)

#FDR with auto corrected P (algorithm GIF)
p1 = ggplot(pv.ppt.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = FDR.sig)) +
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
axis.text.x = element_text(size = 8)
#axis.text.x = element_blank(),
#axis.title.x = element_blank()
)

p2 = ggplot(pv.ppt.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man)) +
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

#FDR correction maunally adjusted P (GIF = 0.85)
ggplot(pv.ppt.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man)) +
#facet_wrap(~scaffold) +
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
#axis.text.x = element_text(angle = 85, size = 10, hjust = 1)
)

pdf("figures/LFMM.ppt.pdf", width = 18, height = 4)
p1
p2
dev.off()


#######################
#Inection duration

#variable for test
X = sample_metadata.site_info$duration_infection

#LFMM ridge
#mod.lfmm = lfmm_ridge(Y = Y, X = X, K = 2)
mod.lfmm.dur_inf = lfmm_ridge(Y = Y, X = X, K = 4)

pv.dur_inf <- lfmm_test(Y = Y,
X = X,
lfmm = mod.lfmm.dur_inf,
calibrate = "gif")

#Example plots
plot(-log10(pv.dur_inf$calibrated.pvalue),
pch = 19,
cex = .3,
xlab = "Probe", ylab = "-Log P",
col = "grey")

plot((pv.dur_inf$B),
pch = 19,
cex = .3,
xlab = "Probe", ylab = "Effect size",
col = "grey")

#Computing genomic inflation factor (GIF) based on calibrated z-scores (http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html) and Francois et al. 2016
lambda = median(pv.dur_inf$score^2)/0.456
lambda #1.09
adj.p.values = pchisq(pv.dur_inf$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
#Note that these calibrated scores are similar as pv$calibrated.pvalue
hist(pv.dur_inf$calibrated.pvalue)

#IN THIS CASE THE CALCULATED VALUES looks conservative

#Try higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.dur_inf$score^2/1.25, df = 1, lower = FALSE)
hist(adj.p.values) #very conservative

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.dur_inf$score^2/1, df = 1, lower = FALSE)
hist(adj.p.values) #This looks good

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv.dur_inf$score^2/0.85, df = 1, lower = FALSE)
hist(adj.p.values)

#THIS IS SHOWING THAT THE GIF CALIBRATION IN THE ALGORITHM IS MORE CONSERVATIVE THAN LOWER VALUES OF LAMBDA
#However, the lower values have a correct distribution under null model
#Try at GIF = 0.1
adj.p.values = pchisq(pv.dur_inf$score^2/1, df = 1, lower = FALSE)
hist(adj.p.values)

#Read scaffold lengths
scf_lens = read.table("data/Nf_post_SPANDx/scaffold_lengths.csv", sep = ",", header = F)
colnames(scf_lens) = c("scaffold", "length")

#Join with actual positiion and chromosome
pv.dur_inf.with_pos = data.frame(calibrated.p = pv.dur_inf$calibrated.pvalue, effect_size = pv.dur_inf$B, SNP_pos, man.adj.p = adj.p.values) %>% left_join(., scf_lens)

#############################
#Outlier tests
# identify the 95% percentile
my_threshold <- quantile((pv.dur_inf.with_pos %>% filter(length > 100000))$calibrated.p, 0.025, na.rm = T)
# make an outlier column in the data.frame
pv.dur_inf.with_pos <- pv.dur_inf.with_pos %>% mutate(outlier = ifelse(calibrated.p < my_threshold, "outlier", "background"))
#Number of outliers
pv.dur_inf.with_pos %>% group_by(outlier) %>% tally()
#Example plot

#3223

#FDR correction
pv.dur_inf.with_pos$FDR.p = p.adjust(pv.dur_inf.with_pos$calibrated.p, method = "fdr", n = length(pv.dur_inf.with_pos$calibrated.p))
pv.dur_inf.with_pos <- pv.dur_inf.with_pos %>% mutate(FDR.sig = ifelse(FDR.p < 0.1, "sig", "background"))
pv.dur_inf.with_pos %>% group_by(FDR.sig) %>% tally()

#23

#FDR correction manual adjustment (lambda = 1)
pv.dur_inf.with_pos$FDR.p.man = p.adjust(pv.dur_inf.with_pos$man.adj.p, method = "fdr", n = length(pv.dur_inf.with_pos$man.adj.p))
pv.dur_inf.with_pos <- pv.dur_inf.with_pos %>% mutate(FDR.sig.man = ifelse(FDR.p.man < 0.05, "sig", "background"))
pv.dur_inf.with_pos %>% group_by(FDR.sig.man) %>% tally()

#23

####################
#ggplots

#Basic plot
ggplot(pv.dur_inf.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p)) +
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
ggplot(pv.dur_inf.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = outlier)) +
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
p1 = ggplot(pv.dur_inf.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = FDR.sig)) +
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

p2 = ggplot(pv.dur_inf.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man)) +
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
