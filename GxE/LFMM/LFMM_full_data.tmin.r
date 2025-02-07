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
source("R_scripts/ggplot_theme.txt")

#The genotype data can simply be read in as a matrix (according the docs)
#OR can try loading LEA and using readLfmm()
Y = as.matrix(read.table("data/Nf_SPANDx_all_seqs/out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode01missing9.lfmm", header = F))

SNP_pos = read.table("data/Nf_SPANDx_all_seqs/out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode01missing9.map")
SNP_pos = SNP_pos[c(1,4)]
colnames(SNP_pos) = c("scaffold", "position")

#principal components analysis for K
pc = prcomp(Y)
plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
points(4,pc$sdev[7]^2, type = "h", lwd = 3, col = "blue")


#read env data. Needs to be sorted by PED .sampleIDs info
source("R_scripts/make_site_metadata.r")

#Join pops to site data

site.info = read.csv("data/sample_metadata/site_info.csv")
site.GDD = read.table("data/sample_metadata/site_climate.GDD.txt", header = T)
site.climate = read.table("data/sample_metadata/sites_climate.txt", header = T)
site.GDD$freezeThaw.annual_mean = site.GDD$freezeThaw.mean_growing + site.GDD$freezeThaw.mean_nongrowing

site_metadata = left_join(site.GDD, site.info %>% select(Site, lat, lon, duration_infection), by = "Site") %>%
    left_join(., site.climate %>% select(Site, tmin, tmax, ppt, MAT, lat, lon, state.name), by = c("Site", "lat", "lon") )


sample_metadata.site_info = left_join(sample_metadata, site_metadata, by = "state.name")
sample_metadata.site_info = sample_metadata.site_info %>% filter(sample %in% fam_info.Nf[,1])
nrow(sample_metadata.site_info)

#######################
#tmin

#variable for test
X = sample_metadata.site_info$tmin

#LFMM ridge
mod.lfmm = lfmm_ridge(Y = Y, X = X, K = 4) #using K = 4 based on PCA and pop structure analyses

pv <- lfmm_test(Y = Y,
X = X,
lfmm = mod.lfmm,
calibrate = "gif")

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
lambda #1.019982
adj.p.values = pchisq(pv$score^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values)
#Note that these calibrated scores are similar as pv$calibrated.pvalue
hist(pv$calibrated.pvalue)
#The auto P values look pretty good

#Try higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv$score^2/1.15, df = 1, lower = FALSE)
hist(adj.p.values)

#Try lower value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv$score^2/0.95, df = 1, lower = FALSE)
hist(adj.p.values)

#Try slightly higher value of GIF -- looking for flat distribution with peak near zero
adj.p.values = pchisq(pv$score^2/1.05, df = 1, lower = FALSE)
hist(adj.p.values)

#THIS IS SHOWING THAT THE GIF CALIBRATION IN THE ALGORITHM IS Good
#However, the (slightly) lower values have a correct distribution under null model
#Try at GIF = 0.95
adj.p.values = pchisq(pv$score^2/0.95, df = 1, lower = FALSE)
hist(adj.p.values)

#Read scaffold lengths
scf_lens = read.table("data/Nf_SPANDx_all_seqs/scaffold_lengths.csv", sep = ",", header = F)
colnames(scf_lens) = c("scaffold", "length")

#Join with actual position and chromosome
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

#5 of 128624 SNPs identified as significant after FDR correction

#FDR correction
#This is based on the manual GIF adjustment
pv.with_pos$FDR.p.man = p.adjust(pv.with_pos$man.adj.p, method = "fdr", n = length(pv.with_pos$man.adj.p))
pv.with_pos <- pv.with_pos %>% mutate(FDR.sig.man = ifelse(FDR.p.man < 0.05, "sig", "background"))
pv.with_pos %>% group_by(FDR.sig.man) %>% tally()

#5 of 128624 SNPs identified as significant after FDR correction

########################
#Rerun the FDR correction to have automatically calculated p-vals
########################
#FDR correction
#This is based on the auto calibartion
pv.with_pos$FDR.p = p.adjust(pv.with_pos$calibrated.p, method = "fdr", n = length(pv.with_pos$calibrated.p))
pv.with_pos <- pv.with_pos %>% mutate(FDR.sig = ifelse(FDR.p < 0.1, "sig", "background"))
pv.with_pos %>% group_by(FDR.sig) %>% tally()

#5 of 130957 SNPs identified as significant after FDR correction

write.table(pv.with_pos, "data/Nf_LFMM_tables/tmin_lfmm.txt", quote = F, row.names = F, sep = "\t")

####################
#ggplots

#Basic plot
ggplot(pv.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p)) +
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
ggplot(pv.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = outlier)) +
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
p1 = ggplot(pv.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = calibrated.p, color = FDR.sig)) +
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
p1

p2 = ggplot(pv.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man)) +
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
p2

#FDR correction maunally adjusted P (GIF = 0.95)
ggplot(pv.with_pos %>% filter(length > 100000), aes(x = position/10^6, y = man.adj.p, color = FDR.sig.man)) +
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

pdf("figures/LFMM.tmin.pdf", width = 18, height = 4)
p1
p2
dev.off()
