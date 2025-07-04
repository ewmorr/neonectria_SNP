library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(scatterpie)
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
        plot.title = element_text(hjust = -0.2, vjust = -2, size = 20) #hjust -0.2 for main, -0.17 for supp
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
        plot.title = element_text(hjust = -0.2, vjust = -2, size = 20) #hjust -0.19 for main, -0.17 for supp
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
    geom_point(position = position_jitter(width = 0.015)) + 
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
        plot.title = element_text(hjust = -0.2, vjust = -2, size = 20) #hjust -0.2 for main, -0.17 for supp
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
    geom_point(position = position_jitter(width = 0.015, height = 0.015)) + 
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
        plot.title = element_text(hjust = -0.205, vjust = -2, size = 20) #hjust -0.19 for main, -0.17 for supp
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
        plot.title = element_text(hjust = -0.20, vjust = -2, size = 20) #hjust -0.19 for main, -0.17 for supp
    )
p5

Nd.pca_scores.metadata %>% filter(Site == "NH.CCM")
####################################
# admixture plots
####################################

#metadata and ancestry 
Nf.fam_info = read.table("data/Nf/final_tables/rm_dups/FINAL_snp.admixture.nosex", header = F)
colnames(Nf.fam_info)[1] = "Sequence_label"
Nf.sample_metadata = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
Nf.sample_metadata = left_join(data.frame(Sequence_label = Nf.fam_info[,1]), Nf.sample_metadata)

Nd.fam_info = read.table("data/Nd/final_tables/rm_dups/FINAL_snp.admixture.nosex", header = F)
colnames(Nd.fam_info)[1] = "Sequence_label"
Nd.sample_metadata = read.csv("data/sample_metadata/Nd_filtered.lat_lon_dur_inf.csv")
Nd.sample_metadata = left_join(data.frame(Sequence_label = Nd.fam_info[,1]), Nd.sample_metadata)


Nf.K_4 = read.table("data/Nf/admixture/FINAL_snp.admixture.4.Q", header = F)
colnames(Nf.K_4) = paste0("Q", 1:ncol(Nf.K_4))
Nf.K_4$Sequence_label = Nf.fam_info$Sequence_label
Nf.K_4$K = 4
Nf.K_4 = Nf.K_4  %>% 
    pivot_longer(cols = starts_with("Q"), names_to = "ancestor", values_to = "Q")

Nd.K_5 = read.table("data/Nd/admixture/FINAL_snp.admixture.5.Q", header = F)
colnames(Nd.K_5) = paste0("Q", 1:ncol(Nd.K_5))
Nd.K_5$Sequence_label = Nd.fam_info$Sequence_label
Nd.K_5$K = 5
Nd.K_5 = Nd.K_5  %>% 
    pivot_longer(cols = starts_with("Q"), names_to = "ancestor", values_to = "Q")


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

Nf.state_order
Nd.state_order

Nf.K_4.clust = left_join(Nf.K_4, Nf_pca_clust) %>%
    left_join(., Nf.pca_scores.metadata %>% select(Sequence_label, MDS1, duration_infection, state) )
Nd.K_5.clust = left_join(Nd.K_5, Nd_pca_clust) %>%
    left_join(., Nd.pca_scores.metadata %>% select(Sequence_label, MDS1, duration_infection, state, Tree_species) )

ancestry_cols = RColorBrewer::brewer.pal(n = 9, "Set1")

##########################################
# structure plots using min centered data
##########################################

#FIRST REDOING ADMIXTURE PLOT WITH NF K=2
Nf.K_2 = read.table("data/Nf/admixture/FINAL_snp.admixture.2.Q", header = F)
colnames(Nf.K_2) = paste0("Q", 1:ncol(Nf.K_2))
Nf.K_2$Sequence_label = Nf.fam_info$Sequence_label
Nf.K_2$K = 2
Nf.K_2 = Nf.K_2  %>% 
    pivot_longer(cols = starts_with("Q"), names_to = "ancestor", values_to = "Q")

Nf.K_2.clust = left_join(Nf.K_2, Nf_pca_clust) %>%
    left_join(., Nf.pca_scores.metadata %>% select(Sequence_label, MDS1, duration_infection, state) )

p6 = ggplot(Nf.K_2.clust, 
                aes(
                    x = reorder(Sequence_label, MDS1), 
                    y = Q, 
                    fill = ancestor
                )
) +
    geom_bar(stat = "identity") +
    facet_grid(
        ~factor(state, 
            levels = Nf.state_order,
        ), 
        scales = "free_x", 
        space = "free_x",
        switch = "x"
    ) + 
    scale_fill_brewer(palette = "Set1", guide = "none") +
    my_gg_theme.def_size +
    labs(
        x = "Site",
        y = "K = 2",
        title = "h"
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        #strip.text.x = element_text(angle = 90, hjust = 1),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = -0.017, vjust = -2.5, size = 20),
        panel.spacing = unit(0.1,'lines'),
        plot.margin = margin(t = -10, r = 5.5, b = -10, l = 5.5)
    )
p6
p7 = ggplot(Nd.K_5.clust, 
                aes(
                    x = reorder(Sequence_label, MDS1), 
                    y = Q, 
                    fill = ancestor
                )
) +
    geom_bar(stat = "identity") +
    facet_grid(
        ~factor(state, 
            levels = Nd.state_order,
        ), 
        scales = "free_x", 
        space = "free_x",
        switch = "x"
    ) + 
    scale_fill_brewer(palette = "Paired", guide = "none") +
    my_gg_theme.def_size +
    labs(
        x = "HCPC cluster (bars represent individuals)",
        y = "K = 5",
        title = "i"
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_blank(),
        #strip.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = -0.055, vjust = -2.5, size = 20),
        panel.spacing = unit(0.1,'lines'),
        plot.margin = margin(t = -10, r = 5.5, b = -10, l = 5.5)
    )
p7



#################################################################
# now moving on with structure plots

#Nf.strK_2 = read.csv("data/Nf/structure/no_locData/K2_centered_on_min.csv")
Nf.strK_2 = read.table("data/Nf/structure/no_locData/K2_ancestry.txt", header = T)

Nf.strK_2$K = 2
Nf.strK_2 = Nf.strK_2  %>% 
    #mutate(across(where(is.numeric), ~ scale(., center = F, scale = sd(.)) )) %>%
    pivot_longer(cols = starts_with("Q"), names_to = "ancestor", values_to = "Q")

#Nd.strK_5 = read.csv("data/Nd/structure/no_locPrior/K5_centered_on_min.csv")
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
                    y = Q, 
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
        title = "j"
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(angle = 90, hjust = 1),
        #strip.text.x = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = -0.015, vjust = -0.25, size = 20),
        panel.spacing = unit(0.1,'lines'),
        plot.margin = margin(t = -10, r = 5.5, b = 5.5, l = 5.5)
    )
p8

nd_str_cols = RColorBrewer::brewer.pal(n = 5, "Paired")

p9 = ggplot(Nd.strK_5.clust, 
                aes(
                    x = reorder(Sequence_label, MDS1), 
                    y = Q, 
                    #y = Q,
                    fill = rev(ancestor)
                )
) +
    geom_bar(stat = "identity", position = "fill") +
    facet_grid(
        ~factor(state, 
            levels = Nd.state_order,
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
        title = "k"
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(angle = 90, hjust = 1),
        #strip.text.x = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = -0.055, vjust = -0.25, size = 20),
        panel.spacing = unit(0.1,'lines'),
        plot.margin = margin(t = -10, r = 5.5, b = 5.5, l = 5.5)
    )
p9


p10 = readRDS("figures/pop_gen/admixture/Nf.avg_map.rds")
p11 = readRDS("figures/pop_gen/admixture/Nd.avg_map.rds")
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
gp11 = ggplotGrob(p11)

gp1 = set_panel_size(g = gp1, width = unit(3, "inches"), height = unit(3, "inches"))
gp2 = set_panel_size(g = gp2, width = unit(3, "inches"), height = unit(3, "inches"))
gp3 = set_panel_size(g = gp3, width = unit(3, "inches"), height = unit(3, "inches"))
gp4 = set_panel_size(g = gp4, width = unit(3, "inches"), height = unit(3, "inches"))
gp5 = set_panel_size(g = gp5, width = unit(3, "inches"), height = unit(3, "inches"))

# we actually don't need to set the panel size for the remaining bc each facet will be treated individually
# across the whole set once they are g-bound
# setting the panel size makes them all *equal* instead of variable
#gp6 = set_panel_size(g = gp6, width = unit(7.7, "inches"))
#gp7 = set_panel_size(g = gp7, width = unit(2.3, "inches"))

gpTop = cbind(gp2, gp1, gp5, gp3, gp4)
#add a few spacer cols to the left of each map
gp10 = gtable_add_cols(gp10, widths = rep(unit(1, "points"),20), pos = 1)
gp11 = gtable_add_cols(gp11, widths = rep(unit(1, "points"),25), pos = 1)

gpTopMiddle = cbind(gp10, gp11)
gpBotMiddle = cbind(gp6, gp7)
gpBottom = cbind(gp8, gp9)

plot(gpTop)
plot(gpTopMiddle)
plot(gpBotMiddle)
plot(gpBottom)

ncol(gpTop)
ncol(gpTopMiddle)
ncol(gpBotMiddle)
ncol(gpBottom)-ncol(gpTop)
#need to add 41 cols
ncol(gpBottom)-ncol(gpTopMiddle)
#need to add 80 cols

gpTop.add = gtable_add_cols(gpTop, widths = rep(unit(0, "points"),20), pos = 1)
ncol(gpTop.add)
gpTop.add = gtable_add_cols(gpTop.add, widths = rep(unit(0, "points"),21), pos = -1)
ncol(gpTop.add)

gpTopMiddle.add = gtable_add_cols(gpTopMiddle, widths = rep(unit(0, "points"),10), pos = 1)
ncol(gpTopMiddle.add)
gpTopMiddle.add = gtable_add_cols(gpTopMiddle.add, widths = rep(unit(0, "points"),25), pos = -1)
ncol(gpTopMiddle.add)


gpAll = rbind(gpTop.add, gpTopMiddle.add, gpBotMiddle, gpBottom)



#pdf("figures/pop_gen/pca.admixture_maps.admixture_structure_site_order.pdf", width = 25.5, height = 14.5)
#grid.arrange(gpTop.add, gpTopMiddle.add, gpBotMiddle, gpBottom, heights = c(0.25,0.4475,0.07,0.11))
#dev.off()


pdf("figures/pop_gen/pca.admixture_maps.admixture_structure_site_order.raw_structure.pdf", width = 25.5, height = 15)
grid.arrange(gpTop.add, gpTopMiddle.add, gpBotMiddle, gpBottom, heights = c(0.25,0.435,0.09,0.125))
dev.off()
