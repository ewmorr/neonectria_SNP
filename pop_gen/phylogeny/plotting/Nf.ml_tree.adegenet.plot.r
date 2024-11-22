library(dplyr)
library(ggplot2)
library(phylogram)
library(ggdendro)
library(dendextend)
library(zoo)
library(pals)
#library(Polychrome)

#############
#Read tree
source("library/ggplot_theme.txt")

fit = readRDS("data/Nf/phylogeny/ML_tree.rds")
fit$tree
Nf.dendro = as.dendrogram.phylo(fit$tree)

#Read metadata
#
sample_metadata = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
colnames(sample_metadata)[4] = "label" #for join below
#
#########


###################################
#pretty plot of tree with ggdendro and using dendextend functions to extract leaves

leaves = Nf.dendro %>% dendextend::get_nodes_attr("leaf")
leaf_height = Nf.dendro %>% dendextend::get_leaves_attr("height")


#######################################
ddata <- ggdendro::dendro_data(Nf.dendro, type = "rectangle")

ddata$labels$leaf_height = leaf_height
ddata$labels = full_join(ddata$labels, sample_metadata)

col_pal = pals::alphabet(n = 18)

p1 = ggplot() +
  geom_segment(data = ddata$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = ddata$labels, aes(x = x, y = leaf_height, color = state, shape = collection_period), size = 4) +
  scale_color_manual(values = c25) +
  scale_shape_manual(values = c(18, 16)) +
  coord_flip() +
 # coord_polar() +
  scale_y_reverse() +
  labs(color = "State/region") +
  theme_dendro() +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.position = c(0.8,0.75)
  )
p1

pdf("figures/pop_gen/phylogeny/Nf_ML_tree.pdf", height = 12, width = 8)
p1
dev.off()


p2 = ggplot() +
  geom_segment(data = ddata$segments, aes(x = x, y = y/10, xend = xend, yend = yend/10)) +
  geom_point(data = ddata$labels, aes(x = x, y = leaf_height/10, color = state, shape = collection_period), size = 3) +
  scale_color_manual(values = c25) +
    scale_shape_manual(values = c(18, 16)) +
  coord_flip() +
  coord_polar() +
  scale_y_reverse() +
  labs(color = "State/region") +
  theme_dendro() +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

pdf("figures/pop_gen/phylogeny/Nf.ML_tree_polar.pdf", height = 10, width = 12)
p2
dev.off()


################################
#ggdendro with coloring of branches
source("library/dendro_data_k.r")

k = 4   # Number of clusters
dendr = dendro_data_k(Nf.dendro, k = k)


#####################
#Add heights and cols to leaves for plotting points
#leaves = Nf.dendro %>% dendextend::get_nodes_attr("leaf")
#leaf_height = Nf.dendro %>% dendextend::get_leaves_attr("height")
#dendr$labels$leaf_height = leaf_height

dendr$labels = full_join(dendr$labels, sample_metadata)


#with points

p3 = ggplot() + 
  geom_segment(data = dendr$segments, 
               aes(x=x, y=y, xend=xend, yend=yend, colour=factor(clust)), 
               lineend = "square", show.legend = FALSE) + 
  scale_colour_manual(values = c("grey60", c25)) +
  geom_point(data = dendr$labels, aes(x = x, y = leaf_height, fill = state, shape = collection_period), size = 4, shape = 21) +
  scale_shape_manual(values = c(23, 21)) +
  scale_fill_manual(values = c25) +
  scale_y_reverse(expand = c(0.2, 0)) + 
  labs(x = NULL, y = NULL, fill = "State/region") +
  coord_flip() +
  theme_dendro() +
  theme(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position = c(0.8, 0.75)
        )
p3

pdf("figures/pop_gen/phylogeny/Nf.ML_tree.k_4.pdf", height = 12, width = 12)
p3
dev.off()




p4 = ggplot() + 
  geom_segment(data = segment(dendr), 
               aes(x=x, y=y, xend=xend, yend=yend, size=factor(line), colour=factor(clust)), 
               lineend = "square", show.legend = FALSE) + 
  scale_colour_manual(values = c("grey60", c25)) +
  scale_size_manual(values = c(.1, 1)) +
  geom_point(data = dendr$labels, aes(x = x, y = leaf_height, fill = state), size = 3, shape = 21) +
  scale_fill_manual(values = c25) +
  scale_y_reverse(expand = c(0.2, 0)) + 
  labs(x = NULL, y = NULL, fill = "State/region") +
  coord_flip() +
  coord_polar() +
  theme_dendro() +
  theme(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)
  )
p4

pdf("figures/Nf.ML_tree_polar.k_4.pdf", height = 10, width = 10)
p4
dev.off()
