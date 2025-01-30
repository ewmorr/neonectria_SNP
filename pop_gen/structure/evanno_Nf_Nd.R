library(ggplot2)
library(dplyr)
library(gridExtra)
source("library/ggplot_theme.txt")

Nf = read.table("data/Nf/structure/no_locData/evanno.r_format.txt", header = T)
Nd = read.table("data/Nd/structure/no_locPrior/evanno.r_format.txt", header = T)

Nf
Nd

p1 = ggplot(Nf %>% filter(K > 1 & K < 10), aes(x = K, y = Delta.K)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = c(2,4,6,8)) +
    labs(title = "a", y = expression(Delta~K)) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.22)
    )
p1

p2 = ggplot(Nd %>% filter(K > 1 & K < 10), aes(x = K, y = Delta.K)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = c(2,4,6,8,10)) +
    labs(title = "b", y = expression(Delta~K)) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.1),
        axis.title.y = element_blank()
    )
p2

pdf("figures/pop_gen/structure/Nf_Nd_evanno.pdf", width = 7, height = 3)
grid.arrange(p1,p2,ncol = 2)
dev.off()
