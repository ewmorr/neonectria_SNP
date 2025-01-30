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
        color = "State/region", 
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
        color = "State/region", 
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

    
# join with cluster data and pca axis 1 for sorting
Nf.K_4.clust = left_join(Nf.K_4, Nf_pca_clust) %>%
    left_join(., Nf.pca_scores.metadata %>% select(Sequence_label, MDS1) )
Nd.K_5.clust = left_join(Nd.K_5, Nd_pca_clust) %>%
    left_join(., Nd.pca_scores.metadata %>% select(Sequence_label, MDS1) )

ancestry_cols = RColorBrewer::brewer.pal(n = 9, "Set1")

p6 = ggplot(Nf.K_4.clust, 
                aes(
                    x = reorder(Sequence_label, MDS1), 
                    y = Q, 
                    fill = ancestor
                )
) +
    geom_bar(stat = "identity") +
    facet_grid(
        ~factor(HCPC.cluster, 
                 levels = c(1,2,3,4), # 4,1,3,2
                 labels = c("1: WV", "2: mixed sites including VA and NH.CW individuals", "3: NH.CW", "4: VA")
        ), 
        scales = "free_x", 
        space = "free_x",
        switch = "x"
    ) + 
    scale_fill_brewer(palette = "Set1", guide = "none") +
    my_gg_theme.def_size +
    labs(
        x = "HCPC cluster (bars represent individuals)",
        y = "K = 4",
        title = "f"
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        #axis.title.y = element_blank(),
        plot.title = element_text(hjust = -0.015, vjust = -0.25)
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
        ~HCPC.cluster, 
        scales = "free_x", 
        space = "free_x",
        switch = "x"
    ) + 
    scale_fill_brewer(palette = "Set2", guide = "none") +
    my_gg_theme.def_size +
    labs(
        x = "HCPC cluster (bars represent individuals)",
        y = "K = 5",
        title = "g"
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        #axis.title.y = element_blank(),
        plot.title = element_text(hjust = -0.06, vjust = -0.25)
    )
p7

gp1 = ggplotGrob(p1)
gp2 = ggplotGrob(p2)
gp3 = ggplotGrob(p3)
gp4 = ggplotGrob(p4)
gp5 = ggplotGrob(p5)
gp6 = ggplotGrob(p6)
gp7 = ggplotGrob(p7)

gp1 = set_panel_size(g = gp1, width = unit(3, "inches"), height = unit(3, "inches"))
gp2 = set_panel_size(g = gp2, width = unit(3, "inches"), height = unit(3, "inches"))
gp3 = set_panel_size(g = gp3, width = unit(3, "inches"), height = unit(3, "inches"))
gp4 = set_panel_size(g = gp4, width = unit(3, "inches"), height = unit(3, "inches"))
gp5 = set_panel_size(g = gp5, width = unit(3, "inches"), height = unit(3, "inches"))

# we actually don't need to set the panel size bc each facet will be treated individually
# across the whole set once they are g-bound
# setting the panel size makes them all *equal* instead of variable
#gp6 = set_panel_size(g = gp6, width = unit(7.7, "inches"))
#gp7 = set_panel_size(g = gp7, width = unit(2.3, "inches"))

gpTop = cbind(gp2, gp1, gp5, gp3, gp4)
gpBottom = cbind(gp6, gp7)

ncol(gpTop)
ncol(gpBottom)
gpBottom.add = gtable_add_cols(gpBottom, widths = rep(unit(1, "points"),5), pos = 1)
ncol(gpBottom.add)
gpBottom.add = gtable_add_cols(gpBottom.add, widths = rep(unit(1, "points"),10), pos = -1)

ncol(gpBottom.add)
ncol(gpTop)

#rbind(gpTop, gbBottom.add)
#gpAll = rbind(gpTop, gbBottom.add)


pdf("figures/pop_gen/pca_and_admixture.pdf", width = 25.5, height = 6.35)
grid.arrange(gpTop, gpBottom.add, heights = c(0.6,0.4))
dev.off()


##########################################
# structure plots using min centered data
##########################################

Nf.K_2 = read.csv("data/Nf/structure/no_locData/K2_centered_on_min.csv")
Nf.K_2 = read.table("data/Nf/structure/no_locData/K2_ancestry.txt", header = T)

Nf.K_2$K = 2
Nf.K_2 = Nf.K_2  %>% 
    #mutate(across(where(is.numeric), ~ scale(., center = F, scale = sd(.)) )) %>%
    pivot_longer(cols = starts_with("Q"), names_to = "ancestor", values_to = "Q")

Nd.K_5 = read.csv("data/Nd/structure/no_locPrior/K5_centered_on_min.csv")
Nd.K_5 = read.table("data/Nd/structure/no_locPrior/K5_ancestry.txt", header = T)
Nd.K_5$K = 5
Nd.K_5 = Nd.K_5  %>% 
    #mutate(across(where(is.numeric), ~ scale(., center = F, scale = sd(.)) )) %>%
    pivot_longer(cols = starts_with("Q"), names_to = "ancestor", values_to = "Q")


# join with cluster data and pca axis 1 for sorting
Nf.K_2.clust = left_join(Nf.K_2, Nf_pca_clust) %>%
    left_join(., Nf.pca_scores.metadata %>% select(Sequence_label, MDS1) )
Nd.K_5.clust = left_join(Nd.K_5, Nd_pca_clust) %>%
    left_join(., Nd.pca_scores.metadata %>% select(Sequence_label, MDS1) )

ancestry_cols = RColorBrewer::brewer.pal(n = 9, "Set1")

p6 = ggplot(Nf.K_2.clust, # %>% filter(ancestor == "Q1"), 
                aes(
                    x = reorder(Sequence_label, MDS1), 
                    y = Q, 
                    fill = ancestor
                )
) +
    geom_bar(stat = "identity") +
    facet_grid(
        ~factor(HCPC.cluster, 
                 levels = c(1,2,3,4), # 4,1,3,2
                 labels = c("1: WV", "2: mixed sites including VA and NH.CW individuals", "3: NH.CW", "4: VA")
        ), 
        scales = "free_x", 
        space = "free_x",
        switch = "x"
    ) + 
    scale_fill_brewer(palette = "Set1", guide = "none") +
    my_gg_theme.def_size +
    labs(
        x = "HCPC cluster (bars represent individuals)",
        y = "K = 4",
        title = "f"
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        #axis.title.y = element_blank(),
        plot.title = element_text(hjust = -0.015, vjust = -0.25)
    )
p6

p7 = ggplot(Nd.K_5.clust, 
                aes(
                    x = reorder(Sequence_label, MDS1), 
                    y = Q, 
                    fill = ancestor
                )
) +
    geom_bar(stat = "identity", position = "fill") +
    facet_grid(
        ~HCPC.cluster, 
        scales = "free_x", 
        space = "free_x",
        switch = "x"
    ) + 
    scale_fill_brewer(palette = "Set2", guide = "none") +
    my_gg_theme.def_size +
    labs(
        x = "HCPC cluster (bars represent individuals)",
        y = "K = 5",
        title = "g"
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        #axis.title.y = element_blank(),
        plot.title = element_text(hjust = -0.06, vjust = -0.25)
    )
p7

