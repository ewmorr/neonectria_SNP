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

sum(1/(seq(from=1, to=(n-1), by=1)))

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
##################################################
##################################################
# the below is not currently able to be used for theta calc
# but instead we can runpixy with a window size of 1
# then we can take the sum of the ratio of count diffs to count comps as the 
# number of seg sites to account for missing data and take the ratio of the sum 
# of count comparisons to number of individuals (i.e., the possible comparisons 
# per window) as the correction factor for sequence length.
# Actually we can probably get to this with the 10Kb window because we have the 
# number of sites with at least one valid gt... can this 2nd ratio be 
# no_sites*number of possible unique combinations (based on sample size)
# 
# No we still need to run this on a per site basis otherwise we don't know if 
# the raw diffs occur at same site or different sites
# then we can use the number of count_comparisons for the harmonic number correction
# but need to solve for sample size
# or can use n minus the count_missing
pi_dat.means = pi_dat %>%
    group_by(pop) %>%
    summarize(
        pi_mean = sum(count_diffs, na.rm = T)/(sum(count_comparisons, na.rm = T)/1000), # as described in pixy docs but per Kb
        diffs_sum = sum(count_diffs, na.rm = T),
        comps_sum = sum(count_comparisons, na.rm = T)
    )

colnames(pi_dat.means)[1] = "state"

pi_dat.meta = left_join(pi_dat.means %>% filter(!state %in% low_n), site_dat) %>%
    left_join(., state_n)

pi_dat.meta$theta_w_per_kb = pi_dat.meta$diffs_sum/harmonic_number(pi_dat.meta$n)/(pi_dat.meta$comps_sum/1000)
pi_dat.meta$theta_w_uncorrected = pi_dat.meta$diffs_sum/harmonic_number(pi_dat.meta$n)

############################
# Some comparisons of different calculations
# 
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
# Actually according to the snpR docs it handles NAs correctly by using
# table of allele frequency at each site (so ignores NAs on a per sample basis).
# It's not clear howvcf-kit is handling NA but it appears it might be replicating
# the vcftools calc (also requires diploid vcf). snpR also provides it's ests
# of theta so is more transparent. Let's go with this for now despite the 
# garbage snpR import functions. (snpR is also published as of 2022)

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



site_cols = read.csv("data/sample_metadata/Nf_site_colors.csv")
site_cols
all_dat.cols = left_join(all_dat, site_cols)


p1 = ggplot(all_dat.cols, aes(x = duration_infection, y = pi_mean)) +
    geom_smooth(method = "lm", color="black") +
    geom_point(aes(fill = state), size = 3, shape = 21) +
    #scale_color_manual(values = all_dat.cols$col) +
    scale_fill_manual(values = c25) +
    labs(
        x = "Duration infestation (years)", 
        y = expression(paste(pi, " (SNPs Kb"^-1, ")"))
    ) +
    my_gg_theme.def_size 

p1 = ggplot(all_dat.cols, aes(x = duration_infection, y = pi_mean)) +
    geom_smooth(method = "lm", color="black") +
    geom_point(aes(fill = state), size = 3, shape = 21) +
    #scale_color_manual(values = all_dat.cols$col) +
    scale_fill_manual(values = all_dat.cols$colors) +
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

p1 = ggplot(all_dat.cols, aes(x = duration_infection, y = pi_mean)) +
    geom_smooth(method = "lm", color="black") +
    geom_point(aes(fill = state), size = 3, shape = 21) +
    #scale_color_manual(values = all_dat.cols$col) +
    scale_fill_manual(values = all_dat.cols$colors, guide = "none") +
    labs(
        x = "Duration infestation (years)", 
        y = expression(paste(pi, " (SNPs Kb"^-1, ")")),
        fill ="Site",
        "a"
    ) +
    annotate(
        geom = "text",
        x = max(all_dat.cols$duration_infection),
        y = min(all_dat.cols$pi_mean),
        hjust = 1,
        label = expression(paste("r = 0.55, ", italic(P), " = 0.04"))
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust= -0.1, vjust = -1)
    )
p1

p2 = ggplot(all_dat.cols, aes(x = duration_infection, y = global_ws.theta/(41018940/1000))) +
    geom_smooth(method = "lm", color="black", linetype = 2) +
    geom_point(aes(fill = state), size = 3, shape = 21) +
    #scale_color_manual(values = all_dat.cols$col) +
    scale_fill_manual(values = all_dat.cols$colors, guide = "none") +
    labs(
        x = "Duration infestation (years)", 
        y = expression(paste(theta[italic(W)], " (SNPs Kb"^-1,")")),
        fill ="Site",
        "b"
    ) +
    annotate(
        geom = "text",
        x = min(all_dat.cols$duration_infection),
        y = max(all_dat.cols$global_ws.theta)/(41018940/1000),
        hjust = 0,
        label = expression(paste("r = 0.18, ", italic(P), " = 0.54"))
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust= -0.1, vjust = -1)
    )
p2

p3 = ggplot(all_dat.cols, aes(x = duration_infection, y = TajimaD.snpR)) +
    geom_smooth(method = "lm", color="black", linetype = 2) +
    geom_point(aes(fill = state), size = 3, shape = 21) +
    #scale_color_manual(values = all_dat.cols$col) +
    scale_fill_manual(values = all_dat.cols$colors) +
    labs(
        x = "Duration infestation (years)", 
        y = "Tajima's D",
        fill ="Site",
        "c"
    ) +
    annotate(
        geom = "text",
        x = min(all_dat.cols$duration_infection),
        y = min(all_dat.cols$TajimaD.snpR),
        hjust = 0,
        label = expression(paste("r = -0.34, ", italic(P), " = 0.23"))
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust= -0.1, vjust = -1)
    )
p3

pdf("figures/pop_gen/pixy/pi_theta_D.snpR_theta_D.pdf", width = 16, height = 4.5)
grid.arrange(p1,p2,p3,widths = c(0.31,0.31,0.38))
dev.off()
