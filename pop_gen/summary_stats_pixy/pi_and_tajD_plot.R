library(dplyr)
library(ggplot2)
library(gridExtra)
source("library/ggplot_theme.txt")

#calc for theta W
harmonic_number = function(x){
    hnv = vector(mode = "numeric", length = length(x))
    for(v in 1:length(x)){
        #hn = 0
        #for(i in 1:x[v]){
        #    hn = hn + 1/i
        #}
        hnv[v] = sum(1/(seq(from=1, to=x[v], by=1)))
    }
    return(hnv)
}
harmonic_number(1:5)

harmonic_number(seq(1:5))


calc_c = function(x){
    c = x*(x-1)/2
    return(c)
}
calc_n = function(x){
    n = (1 + sqrt(1 + 8*x))/2
    return(n)
}
calc_c(7)
calc_n(21)


sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
state_n = sample_metadata.Nf %>% group_by(state) %>% summarise(n = n())
state_n %>% print(n = Inf)
#filter based on n
n_min = 4
low_n = state_n %>% filter(n < n_min) %>% pull(state)

########################
#setting up geo distance
site_dat = sample_metadata.Nf %>% 
    select(state, duration_infection) %>%
    filter(!state %in% low_n) %>% 
    distinct
########################

pi_dat = read.table("data/Nf/pixy/windowed_10kb/pixy_pi.txt", header = T)
theta_dat = read.table("data/Nf/pixy/whole_contig/pixy_watterson_theta.txt", header = T)
tajD_dat = read.table("data/Nf/pixy/whole_contig/pixy_tajima_d.txt", header = T)
##################################################
##################################################
pi_dat.means = pi_dat %>%
    group_by(pop) %>%
    summarize(
        pi_mean = sum(count_diffs, na.rm = T)/(sum(count_comparisons, na.rm = T)/1000), # as described in pixy docs but per Kb
        diffs_sum = sum(count_diffs, na.rm = T),
        comps_sum = sum(count_comparisons, na.rm = T)
    )
theta_dat.means = theta_dat %>%
    group_by(pop) %>%
    summarize(
        theta_mean = sum(raw_watterson_theta, na.rm = T)/(sum(weighted_no_sites, na.rm = T)/1000)
    )
tajD_dat.means = tajD_dat %>%
    group_by(pop) %>%
    summarize(
        tajD_mean = ( sum(raw_pi) - sum(raw_watterson_theta) )/sum(tajima_d_stdev)
    )

colnames(pi_dat.means)[1] = "state"
colnames(theta_dat.means)[1] = "state"
colnames(tajD_dat.means)[1] = "state"

pi_dat.meta = left_join(pi_dat.means %>% filter(!state %in% low_n), theta_dat.means %>% filter(!state %in% low_n)) %>%
    left_join(., tajD_dat.means %>% filter(!state %in% low_n)) %>%
    left_join(., site_dat) %>%
    left_join(., state_n)

#pi_dat.meta$theta_w_per_kb = pi_dat.meta$diffs_sum/harmonic_number(pi_dat.meta$n)/(pi_dat.meta$comps_sum/1000)
#pi_dat.meta$theta_w_uncorrected = pi_dat.meta$diffs_sum/harmonic_number(pi_dat.meta$n)

############################
# Some comparisons of different calculations
# 
# pi
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
# 
# 
# theta
plot(theta_mean ~ duration_infection, pi_dat.meta)
cor(pi_dat.meta$theta_mean, pi_dat.meta$duration_infection)
cor.test(pi_dat.meta$theta_mean, pi_dat.meta$duration_infection)
#r = 0.3388946 
#p-value = 0.2359
#
#
plot(theta_mean ~ duration_infection, pi_dat.meta %>% filter(state != "VA"))
cor(pi_dat.meta$theta_mean[pi_dat.meta$state != "VA"], pi_dat.meta$duration_infection[pi_dat.meta$state != "VA"])
cor.test(pi_dat.meta$theta_mean[pi_dat.meta$state != "VA"], pi_dat.meta$duration_infection[pi_dat.meta$state != "VA"])
#r = -0.01820871  
#p-value = 0.9529
#
#
#tajima's D
plot(tajD_mean ~ duration_infection, pi_dat.meta)
cor(pi_dat.meta$tajD_mean, pi_dat.meta$duration_infection)
cor.test(pi_dat.meta$tajD_mean, pi_dat.meta$duration_infection)
#r = 0.5055843 
#p-value = 0.06513
#
#
plot(tajD_mean ~ duration_infection, pi_dat.meta %>% filter(state != "VA"))
cor(pi_dat.meta$tajD_mean[pi_dat.meta$state != "VA"], pi_dat.meta$duration_infection[pi_dat.meta$state != "VA"])
cor.test(pi_dat.meta$tajD_mean[pi_dat.meta$state != "VA"], pi_dat.meta$duration_infection[pi_dat.meta$state != "VA"])
#r = 0.01949798  
#p-value = 0.9496


tajD.vk = read.csv("data/Nf/pixy/tajD.vcfkit.csv")
colnames(tajD.vk)[2] = "TajimaD.vk"
tajD.snpR = read.csv("data/Nf/pixy/tajD.snpR.csv")
colnames(tajD.snpR)[4] = "TajimaD.snpR"
tajD.snpR
all_dat = left_join(pi_dat.meta, tajD.vk, by = c("state", "duration_infection")) %>%
    left_join(.,tajD.snpR)
#check correaltion between the dif tajD measures
#
#vk versus snpR
plot(TajimaD.vk ~ TajimaD.snpR, all_dat)
cor(all_dat$TajimaD.vk, all_dat$TajimaD.snpR)
cor.test(all_dat$TajimaD.vk, all_dat$TajimaD.snpR)
#pixy vs vk
plot(TajimaD.vk ~ tajD_mean, all_dat)
cor(all_dat$TajimaD.vk, all_dat$tajD_mean)
cor.test(all_dat$TajimaD.vk, all_dat$tajD_mean)
#pixy vs snpR
plot(TajimaD.snpR ~ tajD_mean, all_dat)
cor(all_dat$TajimaD.snpR, all_dat$tajD_mean)
cor.test(all_dat$TajimaD.snpR, all_dat$tajD_mean)
#far from significant correlation in any case
#
# I think we trust pixy because the calculations make sense and are based on 
# established software (skikit-allel and pixy)

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
#r = 0.1782193, P = 0.5421

#snpR ts_theta
plot(global_ts.theta ~ duration_infection, all_dat)
cor(all_dat$global_ts.theta, all_dat$duration_infection)
cor.test(all_dat$global_ts.theta, all_dat$duration_infection)

#snpR n_seg
plot(global_num_seg ~ duration_infection, all_dat)
cor(all_dat$global_num_seg, all_dat$duration_infection)
cor.test(all_dat$global_num_seg, all_dat$duration_infection)

#snpR ws_theta
plot(global_ws.theta ~ all_dat$TajimaD.snpR, all_dat)
cor(all_dat$global_ws.theta, all_dat$TajimaD.snpR)
cor.test(all_dat$global_ws.theta, all_dat$TajimaD.snpR)

#snpR ws_theta
plot(global_ts.theta ~ all_dat$TajimaD.snpR, all_dat)
cor(all_dat$global_ts.theta, all_dat$TajimaD.snpR)
cor.test(all_dat$global_ts.theta, all_dat$TajimaD.snpR)

#note ts.theta is the version ofpi used in the snpR calc
plot(global_ts.theta ~ pi_mean, all_dat)
cor(all_dat$global_ts.theta, all_dat$pi_mean)
cor.test(all_dat$global_ts.theta, all_dat$pi_mean)
#not corellated and pretty ugly


compare_tajD = all_dat %>% 
    select(state, tajD_mean, TajimaD.vk, TajimaD.snpR, n)
colnames(compare_tajD)[2:4] = c("pixy", "vcfKit", "snpR")
compare_tajD$n
plot(compare_tajD[2:4])


site_cols.df = read.csv("data/sample_metadata/Nf_site_colors.csv")
#site_cols
#all_dat.cols = left_join(all_dat, site_cols)

site_cols = site_cols.df$colors
names(site_cols) = site_cols.df$state


p1 = ggplot(all_dat, aes(x = duration_infection, y = pi_mean)) +
    geom_smooth(method = "lm", color="black") +
    geom_point(aes(fill = state), size = 3, shape = 21) +
    #scale_color_manual(values = all_dat.cols$col) +
    scale_fill_manual(values = site_cols) +
    labs(
        x = "Duration infestation (years)", 
        y = expression(paste(pi, " (SNPs Kb"^-1, ")"))
    ) +
    my_gg_theme.def_size 
p1

p1 = ggplot(all_dat, aes(x = duration_infection, y = pi_mean)) +
    geom_smooth(method = "lm", color="black") +
    geom_point(aes(fill = state), size = 3, shape = 21) +
    #scale_color_manual(values = all_dat.cols$col) +
    scale_fill_manual(values = site_cols) +
    labs(
        x = "Duration infestation (years)", 
        y = expression(paste(pi, " (SNPs Kb"^-1, ")")),
        fill ="Site"
    ) +
    annotate(
        geom = "text",
        x = max(all_dat.cols$duration_infection),
        y = min(all_dat.cols$pi_mean),
        hjust = 1,
        label = expression(paste("r = 0.55, ", italic(P), " = 0.04"))
    ) +
    my_gg_theme.def_size 
p1

pdf("figures/pop_gen/pixy/pi_durInf.pdf", width = 6, height = 4.5)
p1
dev.off()

p1 = ggplot(all_dat, aes(x = duration_infection, y = pi_mean)) +
    geom_smooth(method = "lm", color="black") +
    geom_point(aes(fill = state), size = 3, shape = 21) +
    #scale_color_manual(values = all_dat.cols$col) +
    scale_fill_manual(values = site_cols, guide = "none") +
    labs(
        x = "Duration infestation (years)", 
        y = expression(paste(pi, " (SNPs Kb"^-1, ")")),
        fill ="Site",
        "a"
    ) +
    annotate(
        geom = "text",
        x = max(all_dat$duration_infection),
        y = min(all_dat$pi_mean),
        hjust = 1,
        label = expression(paste("r = 0.55, ", italic(P), " = 0.04"))
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust= -0.1, vjust = -1)
    )
p1

p2 = ggplot(all_dat, aes(x = duration_infection, y = theta_mean)) +
    geom_smooth(method = "lm", color="black", linetype = 2) +
    geom_point(aes(fill = state), size = 3, shape = 21) +
    #scale_color_manual(values = all_dat.cols$col) +
    scale_fill_manual(values = site_cols, guide = "none") +
    labs(
        x = "Duration infestation (years)", 
        y = expression(paste(theta[italic(W)], " (SNPs Kb"^-1,")")),
        fill ="Site",
        "b"
    ) +
    annotate(
        geom = "text",
        x = max(all_dat$duration_infection),
        y = min(all_dat$theta_mean),
        hjust = 1,
        label = expression(paste("r = 0.34, ", italic(P), " = 0.24"))
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust= -0.1, vjust = -1)
    )
p2

p3 = ggplot(all_dat, aes(x = duration_infection, y = tajD_mean)) +
    geom_smooth(method = "lm", color="black", linetype = 2) +
    geom_point(aes(fill = state), size = 3, shape = 21) +
    #scale_color_manual(values = all_dat.cols$col) +
    scale_fill_manual(values = site_cols) +
    labs(
        x = "Duration infestation (years)", 
        y = "Tajima's D",
        fill ="Site",
        "c"
    ) +
    annotate(
        geom = "text",
        x = max(all_dat$duration_infection),
        y = min(all_dat$tajD_mean),
        hjust = 1,
        label = expression(paste("r = 0.51, ", italic(P), " = 0.07"))
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust= -0.1, vjust = -1)
    )
p3

pdf("figures/pop_gen/pixy/pi_theta_D.pixy.pdf", width = 14, height = 4.1)
grid.arrange(p1,p2,p3,widths = c(0.315,0.315,0.37))
dev.off()


###############################
###############################
#comparing different tajD calcs
compare_tajD

p1 = ggplot(compare_tajD, aes(x = pixy, y = vcfKit)) +
    geom_point(aes(fill = state), size = 3, shape = 21) +
    scale_fill_manual(values = site_cols, guide = "none") +
    my_gg_theme.def_size 
p1

p2 = ggplot(compare_tajD, aes(x = pixy, y = snpR)) +
    geom_point(aes(fill = state), size = 3, shape = 21) +
    scale_fill_manual(values = site_cols, guide = "none") +
    my_gg_theme.def_size 
p2

p3 = ggplot(compare_tajD, aes(x = snpR, y = vcfKit)) +
    geom_point(aes(fill = state), size = 3, shape = 21) +
    scale_fill_manual(values = site_cols) +
    my_gg_theme.def_size 
p3
    
pdf("figures/pop_gen/pixy/compare_tajD_calcs.pdf", width = 15, height = 4.25)
grid.arrange(p1,p2,p3, widths = c(0.31,0.31,0.38))
dev.off()






