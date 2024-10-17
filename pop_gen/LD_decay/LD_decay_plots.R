library(dplyr)
library(ggplot2)
source("library/ggplot_theme.txt")
library(data.table)
library(RcppRoll)
library(gridExtra)

Nf.r_tab = fread("data/Nf/final_tables/rm_dups/FINAL_snp.1Mb_max.LD_decay_9999.ld.apos-bpos-r", header = F)
Nd.r_tab = fread("data/Nd/final_tables/rm_dups/FINAL_snp.1Mb_max.LD_decay_9999.ld.apos-bpos-r", header = F)
Nc.r_tab = fread("data/Nc/final_tables/FINAL_snp.1Mb_max.LD_decay_9999.ld.apos-bpos-r", header = F)

colnames(Nf.r_tab) = c("posA", "posB", "r2")
colnames(Nd.r_tab) = c("posA", "posB", "r2")
colnames(Nc.r_tab) = c("posA", "posB", "r2")

head(Nf.r_tab)
Nf.r_tab$dist = abs(Nf.r_tab$posA - Nf.r_tab$posB)
Nd.r_tab$dist = abs(Nd.r_tab$posA - Nd.r_tab$posB)
Nc.r_tab$dist = abs(Nc.r_tab$posA - Nc.r_tab$posB)

Nf.r_tab$spp = "Nf"
Nd.r_tab$spp = "Nd"
Nc.r_tab$spp = "Nc"

head(Nf.r_tab)
all.R_tab = rbind(Nf.r_tab, Nd.r_tab, Nc.r_tab)
gc()

#all vals plot (it's an obnoxiously big file)
p1 = ggplot(all.R_tab, aes(x = dist/1000, y = r2)) +
    geom_point(alpha = 0.1) +
    geom_smooth() +
    facet_wrap(~spp, ncol = 3) +
    #scale_x_continuous(limits = c(NA, 100)) +
    labs(x = "Distance (Kb)", y = expression(paste("r"^2))) +
    my_gg_theme

pdf("figures/pop_gen/LD_decay/all_spp.pdf", width = 20, height = 10)
p1
dev.off()

####################
#testing rolling mean v sliding window

Nf.r_tab.sorted = Nf.r_tab[order(Nf.r_tab$dist, decreasing = F),c("r2","dist")]

r2_roll = roll_mean(x = Nf.r_tab.sorted$r2, n = 1000, by = 500)
dist_roll = roll_mean(x = Nf.r_tab.sorted$dist, n = 1000, by = 500)
plot(r2_roll ~ dist_roll)
length(r2_roll)

LD_sliding_window_mean = function(r2,dist_bp,window_size,slide){
    r2.sort = r2[order(dist_bp, decreasing = F)]
    dist.sort = dist_bp[order(dist_bp, decreasing = F)]
    
    i = 1 #initialize start & stop
    i2 = window_size
    
    mean_list = list()
    
    while(i <= (dist.sort[length(dist.sort)] - window_size) + 1){
        mean_list[[i]] = data.frame(
            mean_r2 = mean(r2.sort[dist.sort >= i & dist.sort <= i2], na.rm = T), 
            quant90 = quantile( x = r2.sort[dist.sort >= i & dist.sort <= i2], probs = 0.9, na.rm = T),
            start_index = i,
            stop_index = i2
        )
        i = i+slide
        i2 = i2+slide
    }
    return(dplyr::bind_rows(mean_list))
}

Nf_sliding_window = LD_sliding_window_mean(Nf.r_tab.sorted$r2, Nf.r_tab.sorted$dist, 1000, 500)
plot(Nf_sliding_window$mean_r2 ~ Nf_sliding_window$start_index)

Nf.ecdf = ecdf(Nf.sliding_window$mean_r2)
Nf.ecdf_roll = ecdf(r2_roll)
length(Nf.ecdf)
plot(Nf.ecdf)
plot(Nf.ecdf_roll, add = T, col = "red")

roll_v_window.ks = ks.test(x = Nf_sliding_window$mean_r2, y = r2_roll)
#k.s.test says they are different. We should use the sliding window approach
# to test between spp since it is not dependent on SNP density

Nf_sliding_window = LD_sliding_window_mean(Nf.r_tab.sorted$r2, Nf.r_tab.sorted$dist, 1000, 250)
Nd_sliding_window = LD_sliding_window_mean(Nd.r_tab$r2, Nd.r_tab$dist, 1000, 250)
Nc_sliding_window = LD_sliding_window_mean(Nc.r_tab$r2, Nc.r_tab$dist, 1000, 250)
plot(Nf_sliding_window$mean_r2 ~ Nf_sliding_window$start_index)
plot(Nd_sliding_window$mean_r2 ~ Nd_sliding_window$start_index)
plot(Nc_sliding_window$mean_r2 ~ Nc_sliding_window$start_index)


Nf_roll = data.frame(
    mean_r2 = roll_mean(x = Nf.r_tab$r2[order(Nf.r_tab$dist)], n = 1000, by = 500),
    mean_dist = roll_mean(x = Nf.r_tab$dist[order(Nf.r_tab$dist)], n = 1000, by = 500)
)
Nd_roll = data.frame(
    mean_r2 = roll_mean(Nd.r_tab$r2[order(Nd.r_tab$dist)], n = 1000, by = 500),
    mean_dist = roll_mean(Nd.r_tab$dist[order(Nd.r_tab$dist)], n = 1000, by = 500)
)
Nc_roll = data.frame(
    mean_r2 = roll_mean(Nc.r_tab$r2[order(Nc.r_tab$dist)], n = 1000, by = 500),
    mean_dist = roll_mean(Nc.r_tab$dist[order(Nc.r_tab$dist)], n = 1000, by = 500)
)


#the Nc values get highly variable past 300K (because that's the longest contig presumably
#let's cut everything off at that point
Nf_sliding_window$spp = "Nf"
Nd_sliding_window$spp = "Nd"
Nc_sliding_window$spp = "Nc"

Nf_roll$spp = "Nf"
Nd_roll$spp = "Nd"
Nc_roll$spp = "Nc"

all_slide = rbind(Nf_sliding_window, Nd_sliding_window, Nc_sliding_window)
all_roll = rbind(Nf_roll, Nd_roll, Nc_roll) # retesting the rolling mean to see if it picks up the same
# shape in Nf, it does (duh); stick with sliding window approach
ggplot(all_roll %>% filter(mean_dist < 250000), aes(mean_dist, mean_r2)) +
    geom_point() +
    facet_wrap(~spp, ncol = 1, scales = "free_y") +
    scale_x_continuous(limits = c(NA, 250000))



p1 = ggplot(all_slide %>% filter(start_index < 100000), aes(start_index/1000, mean_r2)) +
    geom_point(alpha = 0.5) +
    facet_wrap(
        ~factor(
            spp, 
            levels = c("Nf", "Nd", "Nc"), 
            labels = c("N. faginata", "N. ditissima", "N. coccinea")
        ), 
        ncol = 1, 
        scales = "free_y",
        strip.position = "right"
    ) +
    geom_text(
        data = data.frame(
            spp = c("Nf", "Nd", "Nc"),
            lab = c("N. faginata", "N. ditissima", "N. coccinea"),
            x = rep(95,3),
            y = c(0.045,0.115, 0.625)
        ),
        aes(x = x, y = y, label = lab)
    ) +
    my_gg_theme.def_size +
    labs(x = "Pairwise SNP distance (Kbp)", y = expression(paste("Linkage disequilibrium (r"^2, ")")), title = "a") +
    theme(
        strip.background = element_blank(),
        strip.placement = "outsdie",
        strip.text = element_blank(),
        plot.title = element_text(hjust = -0.125, vjust = -2)
    ) +
    ggh4x::facetted_pos_scales(
        y = list(
            scale_y_continuous(breaks = seq(0.01,0.05,0.01)),
            scale_y_continuous(breaks = seq(0.04,0.12,0.02)), #this is the only one we change from default
            scale_y_continuous(breaks = seq(0.48,0.64,0.04))
        )
    ) 
p1

range01 <- function(x, ...){(x-min(x))/(max(x)-min(x))}

Nf_sliding_window.lt150 = Nf_sliding_window %>% filter(start_index < 100000)
Nf_sliding_window.lt150$scale_r2 = range01(Nf_sliding_window.lt150$mean_r2)
Nd_sliding_window.lt150 = Nd_sliding_window %>% filter(start_index < 100000)
Nd_sliding_window.lt150$scale_r2 = range01(Nd_sliding_window.lt150$mean_r2)
Nc_sliding_window = Nc_sliding_window[!is.na(Nc_sliding_window$mean_r2) & Nc_sliding_window$mean_r2 > 0.51,] # there are some NA vas at the long end...
Nc_sliding_window.lt150 = Nc_sliding_window %>% filter(start_index < 100000)
Nc_sliding_window.lt150$scale_r2 = range01(Nc_sliding_window.lt150$mean_r2)
all_slide.lt150 = rbind(Nf_sliding_window.lt150, Nd_sliding_window.lt150, Nc_sliding_window.lt150)

ks.test(x = Nf_sliding_window.lt150$scale_r2, y = Nd_sliding_window.lt150$scale_r2)
#D = 0.385, p-value < 2.2e-16
#alternative hypothesis: two-sided
ks.test(x = Nf_sliding_window.lt150$scale_r2, y = Nc_sliding_window.lt150$scale_r2)
#D = 0.17753, p-value = 1.762e-05
#alternative hypothesis: two-sided
ks.test(x = Nd_sliding_window.lt150$scale_r2, y = Nc_sliding_window.lt150$scale_r2)
#D = 0.33063, p-value < 2.2e-16
#alternative hypothesis: two-sided

#negative exponentional fits
Nf.fit <- nls(scale_r2 ~ SSasymp(start_index, Asym, R0, lrc), data = Nf_sliding_window.lt150[,c("scale_r2","start_index")])
summary(Nf.fit)
Nf_sliding_window.lt150$predicted = predict(Nf.fit, Nf_sliding_window.lt150)
ggplot(Nf_sliding_window.lt150, aes(x = start_index, y = scale_r2)) + geom_point() + 
    geom_line(aes(x = start_index, y = predicted), color = "red")

Nd.fit <- nls(scale_r2 ~ SSasymp(start_index, Asym, R0, lrc), data = Nd_sliding_window.lt150[,c("scale_r2","start_index")])
K0 <- min(Nd_sliding_window.lt150$scale_r2)/2
mod0 <- lm(log(scale_r2+0.000000001) ~ start_index, data = Nd_sliding_window.lt150[,c("scale_r2","start_index")])
summary(mod0)
Nd.fit = nls(scale_r2 ~ a * exp(-S * start_index), data = Nd_sliding_window.lt150[,c("scale_r2","start_index")], start = list(a = coef(mod0)[1], S = abs(2*coef(mod0)[2])))
summary(Nd.fit)
Nd_sliding_window.lt150$predicted = predict(Nd.fit, Nd_sliding_window.lt150)
ggplot(Nd_sliding_window.lt150, aes(x = start_index, y = scale_r2)) + geom_point() + 
    geom_line(aes(x = start_index, y = predicted), color = "red")
#The Nd fit sucks with basic nls and doesn't converge with the self-starting alg
##try the roll mean, as there is higher density at the low end
Nd.fit <- nls(mean_r2 ~ SSasymp(mean_dist, Asym, R0, lrc), data = Nd_roll)
mod0 <- lm(log(mean_r2+0.000000001) ~ mean_dist, data = Nd_roll)
Nd.fit = nls(mean_r2 ~ a * exp(-S * mean_dist), data = Nd_roll, start = list(a = coef(mod0)[1], S = abs(2*coef(mod0)[2])))
Nd_roll$predicted = predict(Nd.fit, Nd_roll)
ggplot(Nd_roll, aes(x = mean_dist, y = mean_r2)) + geom_point() + 
    geom_line(aes(x = mean_dist, y = predicted), color = "red")
#even worse

Nc.fit <- nls(scale_r2 ~ SSasymp(start_index, Asym, R0, lrc), data = Nc_sliding_window.lt150[,c("scale_r2","start_index")])
summary(Nc.fit)
Nc_sliding_window.lt150$predicted = predict(Nc.fit, Nc_sliding_window.lt150)
ggplot(Nc_sliding_window.lt150, aes(x = start_index, y = scale_r2)) + geom_point() + 
    geom_line(aes(x = start_index, y = predicted), color = "red")

Nf_half_max_ind = which(Nf_sliding_window.lt150$mean_r2 <= max(Nf_sliding_window.lt150$mean_r2)/2)[1]#we just take the first bc it is  sequential
Nd_half_max_ind = which(Nd_sliding_window.lt150$mean_r2 <= max(Nd_sliding_window.lt150$mean_r2)/2)[1]
Nc_half_max_ind = which(Nc_sliding_window.lt150$mean_r2 <= max(Nc_sliding_window.lt150$mean_r2, na.rm = T)/2)[1]

Nf_half_max_dist = Nf_sliding_window.lt150$start_index[Nf_half_max_ind]
Nd_half_max_dist = Nd_sliding_window.lt150$start_index[Nd_half_max_ind]
Nc_half_max_dist = Nc_sliding_window.lt150$start_index[Nc_half_max_ind]

Nf_half_max_dist
#12251
Nd_half_max_dist
#1751
Nc_half_max_dist
#NA (never reaches it)

p2 = ggplot(all_slide.lt150, 
       aes(
           start_index/1000, 
           scale_r2, 
           color = factor(
               spp, 
               levels = c("Nf", "Nd", "Nc"), 
               labels = c("N. faginata", "N. ditissima", "N. coccinea")
            )
        )
    )+
    geom_point(alpha = 0.2) +
    #stat_smooth(method = "nls", 
    #        formula = y ~ a * exp(-S * x), 
    #        method.args = list(start = list(a = coef(mod0)[1], S = abs(coef(mod0)[2]))),
    #        se = F
    #) +
    geom_smooth(method = loess, method.args = list(degree = 2, span = 0.4))  +
    my_gg_theme.def_size +
    scale_color_brewer(palette = "Set1") +
    labs(x = "Distance (Kbp)", y = "Scaled LD (range 0-1)", title = "b") +
    theme(
        legend.title = element_blank(),
        legend.position = c(0.85, 0.85),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = -0.125, vjust = -2)
    )
p2

pdf("figures/pop_gen/LD_decay/spp_comps.pdf", width = 12, height = 4.5)
grid.arrange(p1,p2,ncol=2)
dev.off()
