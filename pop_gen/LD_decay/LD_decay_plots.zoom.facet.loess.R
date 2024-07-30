require(dplyr)
require(ggplot2)
require(stringr)

source("~/repo/neonectria_genome_reseq_10072020/R_scripts/ggplot_theme.txt")


r_tab_avg <- read.table("~/neonectria_genome_reseq_10072020/LD.1KB_avg.txt", sep = "\t", header = T)

#plots of mean dists
#mean r^2 summarized of 1Kb distance windows


p1 = ggplot(r_tab_avg %>% filter(mid <= 50000), aes(x = mid/10^3, y = mean, color = CHR_A)) +
  geom_line() +
    facet_wrap(~CHR_A, scales = "free") +
  geom_smooth(method = "loess", se = F) +
  scale_color_manual(values = c25) +
    #scale_x_continuous(limits = c(1, 30)) +
  labs(x = "Distance (Kbp)", y = expression(paste("LD r"^2))) +
  my_gg_theme

pdf("LD_decay.1Kb_mean.zoom.facet.loess.pdf", width = 24, height = 16)
p1
dev.off()
