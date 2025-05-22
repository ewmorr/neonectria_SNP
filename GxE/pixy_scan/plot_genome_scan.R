library(dplyr)
library(ggplot2)
source("library/ggplot_theme.txt")

pi_dat = read.table("data/Nf/GxE/pixy/cluster2_noNC_single_pop/pixy_pi.txt", header = T)
theta_dat = read.table("data/Nf/GxE/pixy/cluster2_noNC_single_pop/pixy_watterson_theta.txt", header = T)
tajD_dat = read.table("data/Nf/GxE/pixy/cluster2_noNC_single_pop/pixy_tajima_d.txt", header = T)
rdadapt_select_env = read.csv("data/Nf/GxE/RDA/no_NC.rdadapt_selected_env.csv")
head(rdadapt_select_env)
scf_lens = rdadapt_select_env %>% select(scaffold, length) %>% unique()

head(pi_dat)
pi_dat$pos = floor((pi_dat$window_pos_1 + pi_dat$window_pos_2)/2)
theta_dat$pos = floor((theta_dat$window_pos_1 + theta_dat$window_pos_2)/2)
tajD_dat$pos = floor((tajD_dat$window_pos_1 + tajD_dat$window_pos_2)/2)

pi_dat$metric = "pi"
theta_dat$metric = "theta"
tajD_dat$metric = "tajD"

div_dat = rbind(
    pi_dat %>% select(scaffold = chromosome, pos, metric, value = avg_pi),
    theta_dat %>% select(scaffold = chromosome, pos, metric, value = avg_watterson_theta),
    tajD_dat %>% select(scaffold = chromosome, pos, metric, value = tajima_d)
) %>% left_join(., scf_lens)

p1 = ggplot(
    div_dat %>% filter(length > 10^5), 
    aes(x = pos/10^6, y = value)
) +
    #geom_point(size = 1, alpha = 0.1) +
    geom_line() +
    facet_grid(factor(metric, levels = c("pi", "theta", "tajD"))~scaffold, scales = "free", space = "free_x") +
    scale_x_continuous(breaks = c(0,1,2,3,4,5)) +
    my_gg_theme.def_size +
    theme(
        strip.text.x = element_blank()
    )

#adding SNPs
p2 = ggplot() +
     geom_vline(
        data = rdadapt_select_env %>% filter(length > 10^5), 
        aes(xintercept = position/10^6, color = max_env_predictor),
        linewidth = 0.75, alpha = 0.5
    ) +
    geom_line(
        data = div_dat %>% filter(length > 10^5), 
        aes(x = pos/10^6, y = value)    
    ) +
    scale_color_brewer(palette = "Paired") +   
    facet_grid(factor(metric, levels = c("pi", "theta", "tajD"))~scaffold, scales = "free", space = "free_x") +
    scale_x_continuous(breaks = c(0,1,2,3,4,5)) +
    my_gg_theme.def_size +
    theme(
        strip.text.x = element_blank()
    )


p3 = ggplot() +
     geom_vline(
        data = rdadapt_select_env %>% filter(length > 10^5), 
        aes(xintercept = position/10^6),
        color = "grey", linewidth = 0.75
    ) +
    geom_line(
        data = div_dat %>% filter(length > 10^5), 
        aes(x = pos/10^6, y = value)    
    ) +
    scale_color_brewer(palette = "Paired") +   
    facet_grid(factor(metric, levels = c("pi", "theta", "tajD"))~scaffold, scales = "free", space = "free_x") +
    scale_x_continuous(breaks = c(0,1,2,3,4,5)) +
    my_gg_theme.def_size +
    theme(
        strip.text.x = element_blank()
    )


png("figures/GxE/RDA/no_NC.selected_env_vars.versus_pixy.png", width = 4000, height = 600)
p1
p2
p3
dev.off()
        
        