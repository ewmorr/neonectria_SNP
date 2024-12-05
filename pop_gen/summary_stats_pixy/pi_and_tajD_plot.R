library(dplyr)
library(ggplot2)
source("library/ggplot_theme.txt")

sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
state_n = sample_metadata.Nf %>% group_by(state) %>% summarise(n = n())
state_n %>% print(n = Inf)
#filter based on n
n_min = 1
low_n = state_n %>% filter(n < n_min) %>% pull(state)

########################
#setting up geo distance
site_dat = sample_metadata.Nf %>% 
    select(state, duration_infection) %>%
    filter(!state %in% low_n) %>% 
    distinct
########################

pi_dat = read.table("data/Nf/pixy/windowed_10kb/pixy_pi.txt", header = T)
pi_dat.means = pi_dat %>%
    group_by(pop) %>%
    summarize(pi_mean = sum(count_diffs, na.rm = T)/(sum(count_comparisons, na.rm = T)/1000))
colnames(pi_dat.means)[1] = "state"
pi_dat.meta = left_join(pi_dat.means %>% filter(!state %in% low_n), site_dat)

plot(pi_mean ~ duration_infection, pi_dat.meta)
cor(pi_dat.meta$pi_mean, pi_dat.meta$duration_infection)
cor.test(pi_dat.meta$pi_mean, pi_dat.meta$duration_infection)
#r = 0.5527458 
#p-value = 0.04037
#
#
plot(pi_mean ~ duration_infection, pi_dat.meta %>% filter(state != "VA"))
cor(pi_dat.meta$pi_mean[pi_dat.meta$state != "VA"], pi_dat.meta$duration_infection[pi_dat.meta$state != "VA"])
cor.test(pi_dat.meta$pi_mean[pi_dat.meta$state != "VA"], pi_dat.meta$duration_infection[pi_dat.meta$state != "VA"])
# not sig when remove VA
# r 0.197
# p = 0.5186


tajD.vk = read.csv("data/Nf/pixy/tajD.vcfkit.csv")
colnames(tajD.vk)[2] = "TajimaD.vk"
tajD.snpR = read.csv("data/Nf/pixy/tajD.snpR.csv")
colnames(tajD.snpR)[4] = "TajimaD.snpR"
tajD.snpR
all_dat = left_join(pi_dat.meta, tajD.vk, by = c("state", "duration_infection")) %>%
    left_join(.,tajD.snpR)
#check correaltion between the dif tajD measures
plot(TajimaD.vk ~ TajimaD.snpR, all_dat)
cor(all_dat$TajimaD.vk, all_dat$TajimaD.snpR)
cor.test(all_dat$TajimaD.vk, all_dat$TajimaD.snpR)
#far from significant correlation
# I think we trust vcf-kit bc the TajD is published
# and the snpR package is kind of garbage.
# 

#vk D
plot(TajimaD.vk ~ duration_infection, all_dat)
cor(all_dat$TajimaD.vk, all_dat$duration_infection)
cor.test(all_dat$TajimaD.vk, all_dat$duration_infection)

#snpR D
plot(TajimaD.snpR ~ duration_infection, all_dat)
cor(all_dat$TajimaD.snpR, all_dat$duration_infection)
cor.test(all_dat$TajimaD.snpR, all_dat$duration_infection)

#snpR ws_theta
plot(global_ws.theta ~ duration_infection, all_dat)
cor(all_dat$global_ws.theta, all_dat$duration_infection)
cor.test(all_dat$global_ws.theta, all_dat$duration_infection)

#snpR ts_theta
plot(global_ts.theta ~ duration_infection, all_dat)
cor(all_dat$global_ts.theta, all_dat$duration_infection)
cor.test(all_dat$global_ts.theta, all_dat$duration_infection)

#snpR n_seg
plot(global_num_seg ~ duration_infection, all_dat)
cor(all_dat$global_num_seg, all_dat$duration_infection)
cor.test(all_dat$global_num_seg, all_dat$duration_infection)

site_cols = sample_metadata.Nf %>% 
    select(state, duration_infection) %>%
    distinct
site_cols$col = c25[nrow(site_cols):1]

all_dat.cols = left_join(all_dat, site_cols)

ggplot(all_dat.cols, aes(x = duration_infection, y = pi_mean)) +
    geom_smooth(method = "lm") +
    geom_point(aes(color = state)) +
    #scale_color_manual(values = all_dat.cols$col) +
    scale_color_manual(values = c25) +
    my_gg_theme.def_size 
