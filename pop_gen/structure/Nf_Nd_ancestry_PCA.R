library(vegan)
library(dplyr)
library(ggplot2)
library(gridExtra)
source("library/ggplot_theme.txt")

Nf = read.csv("data/Nf/structure/no_locData/K2_ancestry.csv")
Nd = read.csv("data/Nd/structure/no_locPrior/K5_ancestry.csv")
rownames(Nf) = Nf$Sequence_label
rownames(Nd) = Nd$Sequence_label
Nf = Nf[2:3]
Nd = Nd[2:6]

Nf_metadata = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
Nd_metadata = read.csv("data/sample_metadata/Nd_filtered.lat_lon_dur_inf.csv")

Nf_cols = read.csv("data/sample_metadata/Nf_site_colors.csv")
Nd_cols = read.csv("data/sample_metadata/Nd_site_colors.csv")

Nf_cap = capscale(Nf ~ 1)
Nf_scores = scores(Nf_cap)$sites
Nd_cap = capscale(Nd ~ 1)
Nd_scores = scores(Nd_cap)$sites

# reloading bc pca is kind of useless
Nf = read.csv("data/Nf/structure/no_locData/K2_ancestry.csv")
Nd = read.csv("data/Nd/structure/no_locPrior/K5_ancestry.csv")

Nf_scores.meta = left_join(
    Nf,
    Nf_metadata
)

Nd_scores.meta = left_join(
    Nd,
    Nd_metadata
)


Nf_cols.v = Nf_cols$colors
names(Nf_cols.v) = Nf_cols$state

Nd_cols.v = Nd_cols$colors
names(Nd_cols.v) = Nd_cols$state

p1 = ggplot(Nf_scores.meta, aes(x = Q1, y = Q2, color = state)) +
    geom_point(position = position_jitter(height = 0.0005, width = 0)) +
    scale_color_manual(values = Nf_cols.v) +
    my_gg_theme.def_size
p1    

p2 = ggplot(Nd_scores.meta, aes(x = Q2, y = Q4, color = state)) +
    geom_point(position = position_jitter(height = 0.00025, width = 0.00025)) +
    scale_color_manual(values = Nd_cols.v) +
    my_gg_theme.def_size
p2

#center the ancestry data on the mins
Nf.centered = Nf %>%
    mutate(across(where(is.numeric), ~ .x - min(.)))
    
Nd.centered = Nd %>%
    mutate(across(where(is.numeric), ~ .x - min(.))) 

write.csv(Nf.centered, "data/Nf/structure/no_locData/K2_centered_on_min.csv", row.names = F, quote = F)
write.csv(Nd.centered, "data/Nd/structure/no_locPrior/K5_centered_on_min.csv", row.names = F, quote = F)

Nf_center.meta = left_join(
    Nf.centered,
    Nf_metadata
)

Nd_center.meta = left_join(
    Nd.centered,
    Nd_metadata
)

p1 = ggplot(Nf_center.meta, aes(x = Q1, y = Q2, color = state)) +
    geom_point(position = position_jitter(height = 0.0005, width = 0)) +
    scale_color_manual(values = Nf_cols.v) +
    my_gg_theme.def_size
p1    

    