library(ggtree)
library(ggplot2)
library(RColorBrewer)
source("library/ggplot_theme.txt")

tree = treeio::read.tree("data/Nf/final_tables/rm_dups/FINAL_snp.snps_only.for_phylogeny.ph") 

fit = readRDS("data/Nf/phylogeny/ML_tree.rds")
tree = fit$tree


plot(tree)

sample_metadata = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
Nf_pca_clust = read.csv("data/sample_metadata/Nf_filtered.HCPC_clust.csv")
sample_metadata = sample_metadata[,4:ncol(sample_metadata)]
#colnames(sample_metadata)[1] = "tip.label"

Nf_pca_clust$spp_clust = paste("Nf cluster", Nf_pca_clust$HCPC.cluster)
Nf_pca_clust

which(tree$tip.label %in% c("NG147", "NG148", "NG150", "NG132"))


ggtree(tree)
# the %<+% operator can be used to attach data to a tree plot on the fly
# treeio has implemented a full_join method
ggtree(tree) %<+% sample_metadata

#look at the internal node labels
# note this labels all including tips
ggtree(tree) + geom_text(aes(label=node), hjust= -0.25)
ggtree(tree) + geom_tiplab()

spp_cols = brewer.pal(name = "Set1", n = 3)
cluster_cols = brewer.pal(name = "Paired", n = 9) %>% rev
length(tree$tip.label)

p1 = ggtree(tree) %<+% Nf_pca_clust +
    geom_tippoint(
        aes(color = spp_clust),
        size = 2
    ) +
    #scale_color_manual(values = c(cluster_cols[2:5], cluster_cols[6:9], cluster_cols[1]), na.translate = F) +
    scale_color_manual(values = c25, na.translate = F) +
    labs(color = "State/region") +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    theme(
        legend.text = element_text(size = 10)
    )
p1

p1 = ggtree(tree) %<+% sample_metadata +
    geom_tippoint(
        aes(color = state),
        size = 2
    ) +
    #scale_color_manual(values = c(cluster_cols[2:5], cluster_cols[6:9], cluster_cols[1]), na.translate = F) +
    scale_color_manual(values = c25, na.translate = F) +
    labs(color = "State/region") +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    theme(
        legend.text = element_text(size = 10)
    )
p1


#ggplot2::xlim(0,0.8) +
#ggplot2::ylim(0,144) +
geom_cladelabel(
    node = 180, 
    label = "N. faginata", 
    color = spp_cols[1],
    #angle = 270,
    offset = 0.0102, #how far are the line and label from the tips
    offset.text = 0.01,
    barsize = 1
) +    
    geom_cladelabel(
        node = 147, 
        label = "N. ditissima", 
        color = spp_cols[2],
        #angle = 270,
        offset = 0.016, #how far are the line and label from the tips
        offset.text = 0.01,
        barsize = 1
    ) +    
    geom_cladelabel(
        node = 176, 
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
    
    
pdf("figures/pop_gen/phylogeny/fugr_outgroup.ggtree.pdf", width = 8, height = 9)
p1
dev.off()

