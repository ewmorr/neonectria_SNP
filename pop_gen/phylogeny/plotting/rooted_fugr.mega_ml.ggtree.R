library(ggtree)
library(ggplot2)
library(RColorBrewer)

tree = treeio::read.tree("data/Fugr1_ref/map_against_Fugr_no_core_extract/mega_default_ml.fugr_root.nwk")
plot(tree)


sample_metadata = read.csv("data/sample_metadata/core_fugr.lat_lon_dur_inf.csv")
Nf_pca_clust = read.csv("data/sample_metadata/Nf_filtered.HCPC_clust.csv")
Nd_pca_clust = read.csv("data/sample_metadata/Nd_filtered.HCPC_clust.csv")
colnames(sample_metadata)[1] = "tip.label"

Nf_pca_clust$spp_clust = paste("Nf cluster", Nf_pca_clust$HCPC.cluster)
Nf_pca_clust
Nd_pca_clust$spp_clust = paste("Nd cluster", Nd_pca_clust$HCPC.cluster)
Nd_pca_clust
HCPC_all = rbind(Nf_pca_clust, Nd_pca_clust)

which(tree$tip.label %in% c("NG147", "NG148", "NG150", "NG132"))


ggtree(tree)
# the %<+% operator can be used to attach data to a tree plot on the fly
# treeio has implemented a full_join method
ggtree(tree) %<+% sample_metadata

#look at the internal node labels
# note this labels all including tips
ggtree(tree) + geom_text(aes(label=node), hjust= -0.25)
ggtree(tree) + geom_tiplab()
ggtree(tree) + geom_text2(aes(subset=!isTip, label=node), hjust=1) # get the parent nodes from this

spp_cols = brewer.pal(name = "Set1", n = 3)
cluster_cols = brewer.pal(name = "Paired", n = 9) %>% rev
length(tree$tip.label)

p1 = ggtree(tree) %<+% HCPC_all +
    ggplot2::xlim(0,2) +
    ggplot2::ylim(0,144) +
    geom_cladelabel(
        node = 148, 
        label = "N. faginata", 
        color = spp_cols[1],
        #angle = 270,
        offset = 0.0102, #how far are the line and label from the tips
        offset.text = 0.01,
        barsize = 1
    ) +    
    geom_cladelabel(
        node = 260, 
        label = "N. ditissima", 
        color = spp_cols[2],
        #angle = 270,
        offset = 0.016, #how far are the line and label from the tips
        offset.text = 0.01,
        barsize = 1
    ) +    
    geom_cladelabel(
        node = 256, 
        label = "N. coccinea", 
        color = spp_cols[3],
        #angle = 270,
        offset = 0.0206, #how far are the line and label from the tips
        offset.text = 0.01,
        barsize = 1
    ) +    
    geom_cladelabel(
        node = 144, 
        label = "Fusarium graminearum", 
        color = "black",
        #angle = 270,
        offset = 0, #how far are the line and label from the tips
        #offset.text = 0.01,
        hjust = 0
    ) +
    geom_tippoint(
        aes(
            color = factor(spp_clust, 
                levels = c(paste("Nf cluster", 1:4), paste("Nd cluster", 1:5))
            )
        ), 
        size = 1.4, 
        shape = 15
    ) +
    #scale_color_manual(values = c(cluster_cols[2:5], cluster_cols[6:9], cluster_cols[1]), na.translate = F) +
    scale_color_manual(values = cluster_cols, na.translate = F) +
    labs(color = "HCPC cluster") +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    theme(
        legend.position = c(0.15,0.6),
        legend.text = element_text(size = 10)
    )
p1

pdf("figures/pop_gen/phylogeny/fugr_outgroup.mega_ml.ggtree.pdf", width = 8, height = 9)
p1
dev.off()

