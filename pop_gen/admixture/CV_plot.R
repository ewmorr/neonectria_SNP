library(ggplot2)
library(dplyr)
source("library/ggplot_theme.txt")

Nc = read.table("data/Nc/admixture/CV_by_K.text")
Nf = read.table("data/Nf/admixture/CV_by_K.text")
Nd = read.table("data/Nd/admixture/CV_by_K.text")

Nc.CV = data.frame(
    K = as.integer(sub("\\):", "", sub("\\(K=", "", Nc$V3))),
    CV.error = Nc$V4,
    spp = "Nc"
)

Nd.CV = data.frame(
    K = as.integer(sub("\\):", "", sub("\\(K=", "", Nd$V3))),
    CV.error = Nd$V4,
    spp = "Nd"
)

Nf.CV = data.frame(
    K = as.integer(sub("\\):", "", sub("\\(K=", "", Nf$V3))),
    CV.error = Nf$V4,
    spp = "Nf"
)

all.CV = rbind(Nc.CV, Nd.CV, Nf.CV)

p1 = ggplot(all.CV %>% filter(spp != "Nc"), aes(x = K, y = CV.error)) +
    geom_point() +
    geom_line() +
    facet_wrap(
        ~factor(spp, 
            levels = c("Nf", "Nd"), 
            labels = c("N. faginata", "N. ditissima")
        ), 
        scales = "free_y", 
        nrow = 2
    ) +
    my_gg_theme.def_size +
    labs(y = "Cross validation error")

pdf("figures/pop_gen/admixture/CV_error.pdf", width = 4, height = 3.5)
p1
dev.off()
       
       