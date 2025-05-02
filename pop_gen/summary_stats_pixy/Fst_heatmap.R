library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
source("library/ggplot_theme.txt")

sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
state_n = sample_metadata.Nf %>% group_by(state) %>% summarise(n = n())
state_n %>% print(n = Inf)
#filter based on n
n_min = 4
low_n = state_n %>% filter(n < n_min) %>% pull(state)

#fst.dist = readRDS("data/Nf/pixy/windowed_10kb/fst_dist.window_avg.rds")
fst.dist = readRDS("data/Nf/pixy/whole_contig/fst_dist.rds")
fst.mat = as.matrix(fst.dist) 
fst.mat[upper.tri(fst.mat)] = NA
diag(fst.mat) = NA

fst.mat %>%
    as.data.frame %>% 
    rownames_to_column("pop1") %>%
    pivot_longer(-pop1, names_to = "pop2", values_to = "Fst") -> fst.df
fst.mat

factor_order = rownames(fst.mat)
factor_order2 = c("VA", "NH.CW", "WV", factor_order[!factor_order %in% c("VA", "NH.CW", "WV")])
    
p1 = ggplot(fst.df, # %>% filter(!is.na(Fst)), 
       aes(
           factor(pop1, levels = rev(factor_order)), 
           factor(pop2, levels = factor_order), 
           fill = Fst)
    ) +
    geom_tile() +
    my_gg_theme.def_size +
    scale_fill_gradient(low = "white", high = "red", na.value = "light grey") +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.914, 0.8125),
        legend.background = element_rect(color = "black"),
        axis.text.x = element_text(angle = 55, hjust = 1)
    )
p1

pdf("figures/pop_gen/pixy/Fst_lower_tri.whole_contig_calc.pdf", width = 5, height = 5)
p1
dev.off()

###################################
###################################
fst.dist = readRDS("data/Nf/pixy/whole_contig_clusters_pop/fst_dist.rds")
#fst.dist = readRDS("data/Nf/pixy/whole_contig/fst_dist.rds")
fst.mat = as.matrix(fst.dist) 
fst.mat[upper.tri(fst.mat)] = NA
diag(fst.mat) = NA

fst.mat %>%
    as.data.frame %>% 
    rownames_to_column("pop1") %>%
    pivot_longer(-pop1, names_to = "pop2", values_to = "Fst") -> fst.df

factor_order = rownames(fst.mat)
#factor_order2 = c("VA", "NH.CW", "WV", factor_order[!factor_order %in% c("VA", "NH.CW", "WV")])

fst.df$pop1 = sub("clust_", "cluster ", fst.df$pop1)
fst.df$pop2 = sub("clust_", "cluster ", fst.df$pop2)
factor_order = sub("clust_", "cluster ",factor_order)

p2 = ggplot(fst.df, # %>% filter(!is.na(Fst)), 
            aes(
                factor(pop1, levels = rev(factor_order)), 
                factor(pop2, levels = factor_order), 
                fill = Fst)
) +
    geom_tile() +
    my_gg_theme.def_size +
    scale_fill_gradient(low = "white", high = "red", na.value = "light grey") +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.92, 0.8125),
        legend.background = element_rect(color = "black"),
        axis.text.x = element_text(angle = 55, hjust = 1)
    )
p2

pdf("figures/pop_gen/pixy/Fst_lower_tri.whole_contig_calc.PCA_clusters.pdf", width = 5, height = 5)
p2
dev.off()
