library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(gtable)
library(egg)
source("library/ggplot_theme.txt")



Nf.pca_scores.metadata = read.csv("data/pca/Nf_pca.metadata.csv")
Nd.pca_scores.metadata = read.csv("data/pca/Nd_pca.metadata.csv")
Nf.pca_scores.metadata$clust = paste("Nf cluster", Nf.pca_scores.metadata$clust)
Nd.pca_scores.metadata$clust = paste("Nd cluster", Nd.pca_scores.metadata$clust)

Nf_site_cols = read.csv("data/sample_metadata/Nf_site_colors.csv")
Nd_site_cols = read.csv("data/sample_metadata/Nd_site_colors.csv")
Nf_site_cols.v = Nf_site_cols$colors
names(Nf_site_cols.v) = Nf_site_cols$state
Nd_site_cols.v = Nd_site_cols$colors
names(Nd_site_cols.v) = Nd_site_cols$state

Nf_pca_clust = read.csv("data/sample_metadata/Nf_filtered.HCPC_clust.csv")
Nd_pca_clust = read.csv("data/sample_metadata/Nd_filtered.HCPC_clust.csv")
Nf_pca_clust$spp_clust = paste("Nf cluster", Nf_pca_clust$HCPC.cluster)
Nf_pca_clust
Nd_pca_clust$spp_clust = paste("Nd cluster", Nd_pca_clust$HCPC.cluster)
Nd_pca_clust
HCPC_all = rbind(Nf_pca_clust, Nd_pca_clust)
cluster_cols = RColorBrewer::brewer.pal(name = "Paired", n = 9) 
names(cluster_cols) = levels(as.factor(HCPC_all$spp_clust))


# Nf site
p1 = ggplot(Nf.pca_scores.metadata, 
        aes(MDS1, MDS2, 
            color = state, 
            #shape = factor(collection_period, 
            #    levels = c("early", "modern"), 
            #    labels = c("1960s", "contemporary")
            #)
        )
    ) +
    geom_point(size = 2) + 
    scale_color_manual(values = Nf_site_cols.v) +
    scale_shape_manual(values = c(17,16)) +
    labs(
        color = "Sampling location", 
        shape = "Collection period",
        x = "PCA1 (4.1% variance)",
        y = "PCA2 (3.6% variance)",
        title = "b"
    ) +
    guides(color=guide_legend(ncol=2)) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.19, vjust = -2) #hjust -0.2 for main, -0.17 for supp
    ) 
p1

# Nf cluster
p2 = ggplot(Nf.pca_scores.metadata, 
        aes(MDS1, MDS2, color = as.factor(clust))
    ) +
    geom_point() + 
    scale_color_manual(values = cluster_cols, labels = c(1,2,3,4)) +
    stat_ellipse() +
    #geom_mark_ellipse(expand = 0, aes(fill=clust), show.legend = F) +
    labs(
        color = "Cluster", 
        x = "PCA1 (4.1% variance)",
        y = "PCA2 (3.6% variance)",
        title = "a" #c for main, d for supp, b no scree
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.19, vjust = -2) #hjust -0.19 for main, -0.17 for supp
    )
p2


# Nd site
p3 = ggplot(Nd.pca_scores.metadata, 
        aes(MDS1, MDS2, 
            color = state, 
            #shape = factor(collection_period, 
            #    levels = c("early", "modern"), 
            #    labels = c("1960s", "contemporary")
            #)
        )
    ) +
    geom_point(size = 2, position = position_jitter(width = 0.015)) + 
    scale_color_manual(values = Nd_site_cols.v) +
    scale_shape_manual(values = c(17,16)) +
    labs(
        color = "Sampling location", 
        shape = "Collection period",
        x = paste0("PCA1 (6.5% variance)"),
        y = paste0("PCA2 (4.7% variance)"),
        title = "d"
    ) +
    guides(color=guide_legend(ncol=2)) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.19, vjust = -2) #hjust -0.2 for main, -0.17 for supp
    ) 
p3


# Nd tree sp
p4 = ggplot(Nd.pca_scores.metadata, 
        aes(MDS1, MDS2, 
            color = Tree_species, 
            #shape = factor(collection_period, 
            #    levels = c("early", "modern"), 
            #    labels = c("1960s", "contemporary")
            #)
        )
    ) +
    geom_point(size = 2, position = position_jitter(width = 0.015)) + 
    scale_color_brewer(palette = "Set1") +
    scale_shape_manual(values = c(17,16)) +
    guides(shape = "none") + #comment for supp plot
    labs(
        shape = "Collection period",
        color = "Tree species", 
        x = paste0("PCA1 (6.5% variance)"),
        y = paste0("PCA2 (4.7% variance)"),
        title = "e" #b for main, c for supp
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.19, vjust = -2) #hjust -0.19 for main, -0.17 for supp
    )
p4

# Nd cluster
p5 = ggplot(Nd.pca_scores.metadata, aes(MDS1, MDS2, color = clust)) +
    geom_point() + 
    scale_color_manual(values = cluster_cols, labels = c(1,2,3,4,5)) +
    stat_ellipse() +
    labs(
        color = "Cluster", 
        x = paste0("PCA1 (6.5% variance)"),
        y = paste0("PCA2 (4.7% variance)"),
        title = "c" #c for main, d for supp 
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.19, vjust = -2) #hjust -0.19 for main, -0.17 for supp
    )
p5

####################################
# evanno plots
####################################
Nf.ev = read.table("data/Nf/structure/no_locData/evanno.r_format.txt", header = T)
Nd.ev = read.table("data/Nd/structure/no_locPrior/evanno.r_format.txt", header = T)

p6 = ggplot(Nf.ev %>% filter(K > 1 & K < 10), aes(x = K, y = Delta.K)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = c(2,4,6,8)) +
    labs(title = "f", y = expression(Delta~K)) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.19, vjust = -2)
    )
p6

p7 = ggplot(Nd.ev %>% filter(K > 1 & K < 10), aes(x = K, y = Delta.K)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = c(2,4,6,8,10)) +
    labs(title = "g", y = expression(Delta~K)) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.19, vjust = -2)
    )
p7

##############################
#Structure
##############################
Nf.sample_metadata = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
Nd.sample_metadata = read.csv("data/sample_metadata/Nd_filtered.lat_lon_dur_inf.csv")

################
# sorting
# join with cluster data and pca axis 1 for sorting
#sort by infection duration then by state
Nf.try_order = Nf.sample_metadata[with(Nf.pca_scores.metadata, order(duration_infection, state)),"Sequence_label"] 
Nd.try_order = Nd.sample_metadata[with(Nd.pca_scores.metadata, order(duration_infection, state)),"Sequence_label"] 

# ordering states
Nf_states_dur = Nf.sample_metadata %>%
    select(state, duration_infection, Year_infection_observed.min, lat, lon, collection_period, county) %>%
    distinct
Nd_states_dur = Nd.sample_metadata %>%
    select(state, duration_infection, Year_infection_observed.min, lat, lon, collection_period, county) %>%
    distinct

#pushing WV and NC to plot next to PA.W and QC sites will then plot with MI. 
# this makes more sense geographically and the infest age dif is only 2 years
# 
Nf_states_dur[Nf_states_dur$state == "WV" | Nf_states_dur$state == "NC", "Year_infection_observed.min"] = 1996
# for Nd we push QC.OUG to accompliush the same
Nd_states_dur[Nd_states_dur$state == "QC.OUG", "Year_infection_observed.min"] = 2000

Nf.state_order = Nf_states_dur$state[
    order(Nf_states_dur$Year_infection_observed.min, -Nf_states_dur$lon)
]
Nd.state_order = Nd_states_dur$state[
    order(Nd_states_dur$Year_infection_observed.min, -Nd_states_dur$lon)
]

# ancestry data
Nf.strK_2 = read.table("data/Nf/structure/no_locData/K2_ancestry.txt", header = T)

Nf.strK_2$K = 2
Nf.strK_2 = Nf.strK_2  %>% 
    #mutate(across(where(is.numeric), ~ scale(., center = F, scale = sd(.)) )) %>%
    pivot_longer(cols = starts_with("Q"), names_to = "ancestor", values_to = "Q")

Nd.strK_5 = read.table("data/Nd/structure/no_locPrior/K5_ancestry.txt", header = T)
Nd.strK_5$K = 5
Nd.strK_5 = Nd.strK_5  %>% 
    #mutate(across(where(is.numeric), ~ scale(., center = F, scale = sd(.)) )) %>%
    pivot_longer(cols = starts_with("Q"), names_to = "ancestor", values_to = "Q")


# join with cluster data and pca axis 1 for sorting
Nf.strK_2.clust = left_join(Nf.strK_2, Nf_pca_clust) %>%
    left_join(., Nf.pca_scores.metadata %>% select(Sequence_label, MDS1, state, duration_infection) )
Nd.strK_5.clust = left_join(Nd.strK_5, Nd_pca_clust) %>%
    left_join(., Nd.pca_scores.metadata %>% select(Sequence_label, MDS1, state, duration_infection, Tree_species) )

ancestry_cols = RColorBrewer::brewer.pal(n = 3, "Set1")

p8 = ggplot(Nf.strK_2.clust, # %>% filter(ancestor == "Q1"), 
                aes(
                    x = reorder(Sequence_label, MDS1), 
                    y = log(Q*100+1), 
                    #y = Q,
                    fill = ancestor
                )
) +
    geom_bar(stat = "identity",position = "fill") +
    facet_grid(
        ~factor(state, 
            levels = Nf.state_order,
        ), 
        scales = "free_x", 
        space = "free_x",
        switch = "x"
    ) + 
    scale_fill_manual(values = ancestry_cols[2:1], guide = "none") +
    my_gg_theme.def_size +
    labs(
        y = "K = 2",
        title = "h"
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(angle = 90, hjust = 1),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = -0.015, vjust = -0.25),
        panel.spacing = unit(0.0,'lines')
    )
p8

nd_str_cols = RColorBrewer::brewer.pal(n = 5, "Set2")

p9 = ggplot(Nd.strK_5.clust, 
                aes(
                    x = reorder(Sequence_label, MDS1), 
                    y = log(Q*100+1), 
                    #y = Q,
                    fill = rev(ancestor)
                )
) +
    geom_bar(stat = "identity", position = "fill") +
    facet_grid(
        ~factor(state, 
            levels = Nd.state_order[c(1:4,16,5,7,6,8:10,17,18,12,11,14,13,15)]
        ), 
        scales = "free_x", 
        space = "free_x",
        switch = "x"
    ) + 
    scale_fill_manual(values = nd_str_cols[c(5,4,3,1,2)], guide = "none") +
    #scale_fill_brewer(palette = "Set2", guide = "none") +
    my_gg_theme.def_size +
    labs(
        y = "K = 5",
        title = "i"
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(angle = 90, hjust = 1),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = -0.05, vjust = -0.25),
        panel.spacing = unit(0.0,'lines')
    )
p9

p10 = ggplot(Nd.strK_5.clust, 
                aes(
                    x = reorder(Sequence_label, MDS1), 
                    y = log(Q*100+1), 
                    #y = Q,
                    fill = rev(ancestor)
                )
) +
    geom_bar(stat = "identity", position = "fill") +
    facet_grid(
        ~factor(Tree_species
        ), 
        scales = "free_x", 
        space = "free_x",
        switch = "x"
    ) + 
    scale_fill_manual(values = nd_str_cols[c(5,4,3,1,2)], guide = "none") +
    #scale_fill_brewer(palette = "Set2", guide = "none") +
    my_gg_theme.def_size +
    labs(
        y = "K = 5",
        title = "j"
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(angle = 90, hjust = 1),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = -0.05, vjust = -0.25),
        panel.spacing = unit(0.0,'lines')
    )
p10


# setting up big plot
gp1 = ggplotGrob(p1)
gp2 = ggplotGrob(p2)
gp3 = ggplotGrob(p3)
gp4 = ggplotGrob(p4)
gp5 = ggplotGrob(p5)
gp6 = ggplotGrob(p6)
gp7 = ggplotGrob(p7)
gp8 = ggplotGrob(p8)
gp9 = ggplotGrob(p9)
gp10 = ggplotGrob(p10)

gp1 = set_panel_size(g = gp1, width = unit(3, "inches"), height = unit(3, "inches"))
gp2 = set_panel_size(g = gp2, width = unit(3, "inches"), height = unit(3, "inches"))
gp3 = set_panel_size(g = gp3, width = unit(3, "inches"), height = unit(3, "inches"))
gp4 = set_panel_size(g = gp4, width = unit(3, "inches"), height = unit(3, "inches"))
gp5 = set_panel_size(g = gp5, width = unit(3, "inches"), height = unit(3, "inches"))
gp6 = set_panel_size(g = gp6, width = unit(3, "inches"), height = unit(3, "inches"))
gp7 = set_panel_size(g = gp7, width = unit(3, "inches"), height = unit(3, "inches"))

# we actually don't need to set the panel size bc each facet will be treated individually
# across the whole set once they are g-bound
# setting the panel size makes them all *equal* instead of variable
#gp6 = set_panel_size(g = gp6, width = unit(7.7, "inches"))
#gp7 = set_panel_size(g = gp7, width = unit(2.3, "inches"))

gpTop = cbind(gp2, gp1, gp5, gp3, gp4, gp6, gp7)
gpBottom = cbind(gp8, gp9, gp10)

ncol(gpTop)
ncol(gpBottom)
ncol(gpBottom)-ncol(gpTop)

gpTop.add = gtable_add_cols(gpTop, widths = rep(unit(1, "points"),22), pos = 1)
ncol(gpTop.add)
gpTop.add = gtable_add_cols(gpTop.add, widths = rep(unit(1, "points"),22), pos = -1)


ncol(gpTop.add)
ncol(gpBottom)

#rbind(gpTop, gbBottom.add)
#gpAll = rbind(gpTop, gbBottom.add)


pdf("figures/pop_gen/pca_structure_evanno.site_order.pdf", width = 33, height = 8)
grid.arrange(gpTop.add, gpBottom, heights = c(0.6,0.4))
dev.off()
