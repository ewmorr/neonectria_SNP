require(dplyr)
require(ggplot2)
require(stringr)

source("~/repo/neonectria_genome_reseq_10072020/R_scripts/ggplot_theme.txt")

r_tab = read.table("~/Nf_SPANDx_all_seqs/Outputs/Master_vcf/out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.LD_decay.ld.gz", header = T)

#remove uneeded cols to save mem
r_tab = subset(r_tab, select = -c(SNP_A, CHR_B, SNP_B))
#calculate SNP dist
r_tab$distance = r_tab$BP_B - r_tab$BP_A

#plot distance by all points with smoothing
#all points not summarized by dist

#p1 = ggplot(r_tab, aes(x = distance/10^3, y = R2)) +
#  geom_point(alpha = 0.15) +
#  geom_smooth(method = "loess") +
#  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
#  my_gg_theme

#all points not summarized by dist
#colored by chromosome

#p2 = ggplot(r_tab, aes(x = distance/10^3, y = R2, color = CHR_A)) +
#  geom_point(alpha = 0.15) +
#  geom_smooth(method = "loess") +
#  scale_color_manual(values = c25) +
#  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
#  my_gg_theme

#pdf("LD_decay.all_points.pdf", width = 16, height = 6)
#p1
#p2
#dev.off()

#rm(list = c(p1, p2))

###############
#calculate means
#different dists at 1K, 2.7.5K, 7.5K, 10K
###############

#############
#1KB mean
#############
print("starting 1K means")
summarize_dist = 1000
r_tab$distc <- cut(
  r_tab$distance,
  breaks=seq(
    from=0,
    to=ceiling(max(r_tab$distance)/summarize_dist)*summarize_dist,
    by=summarize_dist
    ) 
  )

r_tab_avg <- r_tab %>% group_by(distc, CHR_A) %>% summarise(mean=mean(R2),median=median(R2))

#calc mean distance
r_tab_avg <- r_tab_avg %>% mutate(
    start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
    end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
    mid=start+((end-start)/2)
    )

#plots of mean dists
#mean r^2 summarized of 1Kb distance windows
p3 = ggplot(r_tab_avg, aes(x = mid/10^3, y = mean)) +
  geom_point() +
  geom_smooth(method = "loess", se = F) +
  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
  my_gg_theme

#mean r^2 summarized of 1Kb distance windows
#colored by chromosome
p4 = ggplot(r_tab_avg, aes(x = mid/10^3, y = mean, color = CHR_A)) +
  geom_line() +
  #geom_smooth(method = "loess", se = F) +
  scale_color_manual(values = c25) +
  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
  my_gg_theme

#all points not summarized by dist
#line of mean r^2 summarized of 1Kb distance windows
p5 = ggplot() +
  #geom_point(data = r_tab, aes(x = distance/10^3, y = R2)) +
  geom_line(data = r_tab_avg, aes(x = mid/10^3, y = mean) ) +
  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
  my_gg_theme

pdf("LD_decay.1Kb_mean.pdf", width = 16, height = 6)
p3
p4
p5
dev.off()

rm(r_tab_avg)

#############
#2.5KB mean
#############
print("starting 2.5K means")
summarize_dist = 2500
r_tab$distc <- cut(
  r_tab$distance,
  breaks=seq(
    from=0,
    to=ceiling(max(r_tab$distance)/summarize_dist)*summarize_dist,
    by=summarize_dist
  ) 
)

r_tab_avg <- r_tab %>% group_by(distc, CHR_A) %>% summarise(mean=mean(R2),median=median(R2))

#calc mean distance
r_tab_avg <- r_tab_avg %>% mutate(
  start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
  end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
  mid=start+((end-start)/2)
)

#plots of mean dists
#mean r^2 summarized of 2.5Kb distance windows
p3 = ggplot(r_tab_avg, aes(x = mid/10^3, y = mean)) +
  geom_point() +
  geom_smooth(method = "loess", se = F) +
  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
  my_gg_theme

#mean r^2 summarized of 2.5Kb distance windows
#colored by chromosome
p4 = ggplot(r_tab_avg, aes(x = mid/10^3, y = mean, color = CHR_A)) +
  geom_line() +
  #geom_smooth(method = "loess") +
  scale_color_manual(values = c25) +
  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
  my_gg_theme

#all points not summarized by dist
#line of mean r^2 summarized of 2.5Kb distance windows
p5 = ggplot() +
  #geom_point(data = r_tab, aes(x = distance/10^3, y = R2)) +
  geom_line(data = r_tab_avg, aes(x = mid/10^3, y = mean), color = "black" ) +
  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
  my_gg_theme

pdf("LD_decay.2.5Kb_mean.pdf", width = 16, height = 6)
p3
p4
p5
dev.off()

rm(r_tab_avg)

#############
#5KB mean
#############
print("starting 5K means")
summarize_dist = 5000
r_tab$distc <- cut(
  r_tab$distance,
  breaks=seq(
    from=0,
    to=ceiling(max(r_tab$distance)/summarize_dist)*summarize_dist,
    by=summarize_dist
  ) 
)

r_tab_avg <- r_tab %>% group_by(distc, CHR_A) %>% summarise(mean=mean(R2),median=median(R2))

#calc mean distance
r_tab_avg <- r_tab_avg %>% mutate(
  start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
  end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
  mid=start+((end-start)/2)
)

#plots of mean dists
#mean r^2 summarized of 5Kb distance windows
p3 = ggplot(r_tab_avg, aes(x = mid/10^3, y = mean)) +
  geom_point() +
  geom_smooth(method = "loess", se = F) +
  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
  my_gg_theme

#mean r^2 summarized of 5Kb distance windows
#colored by chromosome
p4 = ggplot(r_tab_avg, aes(x = mid/10^3, y = mean, color = CHR_A)) +
  geom_line() +
  #geom_smooth(method = "loess") +
  scale_color_manual(values = c25) +
  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
  my_gg_theme

#all points not summarized by dist
#line of mean r^2 summarized of 7.5Kb distance windows
p5 = ggplot() +
  #geom_point(data = r_tab, aes(x = distance/10^3, y = R2)) +
  geom_line(data = r_tab_avg, aes(x = mid/10^3, y = mean), color = "black" ) +
  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
  my_gg_theme

pdf("LD_decay.5Kb_mean.pdf", width = 16, height = 6)
p3
p4
p5
dev.off()

rm(r_tab_avg)

#############
#7.5KB mean
#############
print("starting 7.5K means")
summarize_dist = 7500
r_tab$distc <- cut(
  r_tab$distance,
  breaks=seq(
    from=0,
    to=ceiling(max(r_tab$distance)/summarize_dist)*summarize_dist,
    by=summarize_dist
  ) 
)

r_tab_avg <- r_tab %>% group_by(distc, CHR_A) %>% summarise(mean=mean(R2),median=median(R2))

#calc mean distance
r_tab_avg <- r_tab_avg %>% mutate(
  start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
  end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
  mid=start+((end-start)/2)
)

#plots of mean dists
#mean r^2 summarized of 7.5Kb distance windows
p3 = ggplot(r_tab_avg, aes(x = mid/10^3, y = mean)) +
  geom_point() +
  geom_smooth(method = "loess", se = F) +
  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
  my_gg_theme

#mean r^2 summarized of 7.5Kb distance windows
#colored by chromosome
p4 = ggplot(r_tab_avg, aes(x = mid/10^3, y = mean, color = CHR_A)) +
  geom_line() +
  #geom_smooth(method = "loess") +
  scale_color_manual(values = c25) +
  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
  my_gg_theme

#all points not summarized by dist
#line of mean r^2 summarized of 7.5Kb distance windows
p5 = ggplot() +
  #geom_point(data = r_tab, aes(x = distance/10^3, y = R2)) +
  geom_line(data = r_tab_avg, aes(x = mid/10^3, y = mean), color = "black" ) +
  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
  my_gg_theme

pdf("LD_decay.7.5Kb_mean.pdf", width = 16, height = 6)
p3
p4
p5
dev.off()

rm(r_tab_avg)

#############
#10KB mean
#############
print("starting 10K means")
summarize_dist = 10000
r_tab$distc <- cut(
  r_tab$distance,
  breaks=seq(
    from=0,
    to=ceiling(max(r_tab$distance)/summarize_dist)*summarize_dist,
    by=summarize_dist
  ) 
)

r_tab_avg <- r_tab %>% group_by(distc, CHR_A) %>% summarise(mean=mean(R2),median=median(R2))

#calc mean distance
r_tab_avg <- r_tab_avg %>% mutate(
  start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
  end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
  mid=start+((end-start)/2)
)

#plots of mean dists
#mean r^2 summarized of 10Kb distance windows
p3 = ggplot(r_tab_avg, aes(x = mid/10^3, y = mean)) +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
  my_gg_theme

#mean r^2 summarized of 10Kb distance windows
#colored by chromosome
p4 = ggplot(r_tab_avg, aes(x = mid/10^3, y = mean, color = CHR_A)) +
  geom_line() +
  #geom_smooth(method = "loess") +
  scale_color_manual(values = c25) +
  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
  my_gg_theme

#all points not summarized by dist
#line of mean r^2 summarized of 10Kb distance windows
p5 = ggplot() +
  #geom_point(data = r_tab, aes(x = distance/10^3, y = R2)) +
  geom_line(data = r_tab_avg, aes(x = mid/10^3, y = mean), color = "black" ) +
  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
  my_gg_theme

pdf("LD_decay.10Kb_mean.pdf", width = 16, height = 6)
p3
p4
p5
dev.off()
