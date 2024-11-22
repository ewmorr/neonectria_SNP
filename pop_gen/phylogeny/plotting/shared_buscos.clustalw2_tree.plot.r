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

fit = treeio::read.tree("data/shared_buscos/final_tables/rm_dups/FINAL_snp.snps_only.for_phylogeny.ph")
fit
dendro = as.dendrogram.phylo(fit)
plot(dendro)

#Read metadata
#
sample_metadata = read.csv("data/sample_metadata/shared_buscos.lat_lon_dur_inf.csv")
colnames(sample_metadata)
colnames(sample_metadata)[1] = "label" #for join below
#
#########

# first plot

dend.labels.new = left_join(data.frame(label = dendro %>% labels), sample_metadata)
dend.labels.new$label.new = paste0(dend.labels.new$spp, 1:nrow(dend.labels.new))
dend.labels.new$color = vector(mode = "character", length = nrow(dend.labels.new))
dend.labels.new[dend.labels.new$spp == "Nf", "color"] = "blue"
dend.labels.new[dend.labels.new$spp == "Nc", "color"] = "red"
dend.labels.new[dend.labels.new$spp == "Nd", "color"] = "green"

dendro %>% set("labels", dend.labels.new$label.new) %>% plot()

#"raising" doesn't do much for us
dendro %>% 
    set("labels", dend.labels.new$label.new) %>%
    set("labels_colors", dend.labels.new$color) %>%
    raise.dendrogram(0) %>%
    plot()

#pruning
Nd.labels = dend.labels.new %>% filter(spp == "Nd") %>% pull(label.new)
dendro %>% 
    set("labels", dend.labels.new$label.new) %>%
    set("labels_colors", dend.labels.new$color) %>%
    prune(Nd.labels) %>%
    #ladderize() %>%
    plot(horiz = T)

dendro.no_Nd = dendro %>% 
    set("labels", dend.labels.new$label.new) %>%
    prune(Nd.labels) 
dendro.no_Nd %>% plot

dendro %>% 
    set("labels", dend.labels.new$label.new) %>%
    set("labels_colors", dend.labels.new$color) %>%
    #prune(Nd.labels) %>%
    rank_branches %>%
    dendextend::ladderize %>%
    plot(horiz = T)

dendro %>% 
    set("labels", dend.labels.new$label.new) %>%
    set("labels_colors", dend.labels.new$color) %>%
    #prune(Nd.labels) %>%
    #rank_branches %>%
    ladderize %>%
    plot(horiz = T)

dendro %>% 
    set("labels", dend.labels.new$label.new) %>%
    set("labels_colors", dend.labels.new$color) %>%
    prune(Nd.labels) %>%
    #rank_branches %>%
    ladderize %>%
    circlize_dendrogram() %>%
    plot(horiz = T)

dendro %>% 
    set("labels", dend.labels.new$label.new) %>%
    set("labels_colors", dend.labels.new$color) %>%
    prune(Nd.labels) %>%
    rank_branches %>%
    ladderize %>%
    circlize_dendrogram()

dendro %>% 
    set("labels", dend.labels.new$label.new) %>%
    set("labels_colors", dend.labels.new$color) %>%
    #prune(Nd.labels) %>%
    rank_branches %>%
    ladderize %>%
    circlize_dendrogram()

dendro %>% 
    set("labels", dend.labels.new$label.new) %>%
    set("labels_colors", dend.labels.new$color) %>%
#    prune(Nd.labels) %>%
    rank_branches %>%
    ladderize %>%
#    hang.dendrogram() %>%
    circlize_dendrogram()

dendro %>% 
    set("labels", dend.labels.new$label.new) %>%
    set("leaves_pch", 19) %>%
    set("leaves_col", dend.labels.new$color) %>%
    #    prune(Nd.labels) %>%
    rank_branches %>%
    ladderize %>%
#    hang.dendrogram() %>%
    circlize_dendrogram(labels = F)

dendro %>% 
    set("labels", dend.labels.new$label.new) %>%
    set("labels_colors", dend.labels.new$color) %>%
    #set("leaves_pch", 19) %>%
    #set("leaves_col", dend.labels.new$color) %>%
    #prune(Nd.labels) %>%
    rank_branches %>%
    ladderize %>%
#    hang.dendrogram() %>%
    plot(horiz = T)


#try rotation
Nf.labs = dend.labels.new %>% filter(spp == "Nf") %>% pull(label.new)
Nd.labs = dend.labels.new %>% filter(spp == "Nd") %>% pull(label.new)
Nc.labs = dend.labels.new %>% filter(spp == "Nc") %>% pull(label.new)

new_lab_order = c(Nc.labs, Nf.labs, Nd.labs)

dendro %>% 
    set("labels", dend.labels.new$label.new) %>%
    set("labels_colors", dend.labels.new$color) %>%
    #set("leaves_pch", 19) %>%
    #set("leaves_col", dend.labels.new$color) %>%
    prune(Nd.labels) %>%
    rotate(new_lab_order) %>%
    rank_branches %>%
    ladderize %>%
    #    hang.dendrogram() %>%
    plot(horiz = T)

###################################
#pretty plot of tree with ggdendro and using dendextend functions to extract leaves

leaves = dendro.no_Nd %>% dendextend::get_nodes_attr("leaf")
leaf_height = dendro.no_Nd %>% dendextend::get_leaves_attr("height")


#######################################
#dendro = as.dendrogram.phylo(fit)
ddata <- ggdendro::dendro_data(dendro.no_Nd, type = "rectangle")

ddata$labels$leaf_height = leaf_height
colnames(dend.labels.new)[c(1,8)] = c("label.old", "label")
ddata$labels = full_join(ddata$labels, dend.labels.new)

#col_pal = pals::alphabet(n = 18)

p1 = ggplot() +
  geom_segment(data = ddata$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = ddata$labels, aes(x = x, y = leaf_height, color = spp), size = 1) +
  #scale_color_manual(values = c25) +
  scale_shape_manual(values = c(18, 16)) +
  coord_flip() +
 # coord_polar() +
  scale_y_reverse() +
  labs(color = "Species") +
  theme_dendro() +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.position = c(0.8,0.75)
  )
p1

pdf("figures/pop_gen/phylogeny/shared_buscos_ML_tree.pdf", height = 12, width = 8)
p1
dev.off()


p2 = ggplot() +
  geom_segment(data = ddata$segments, aes(x = x, y = y/100, xend = xend, yend = yend/100)) +
  geom_point(data = ddata$labels, aes(x = x, y = leaf_height/100, color = spp, shape = collection_period), size = 1) +
  #scale_color_manual(values = c25) +
    scale_shape_manual(values = c(18, 16)) +
  coord_flip() +
  coord_polar() +
  scale_y_reverse() +
  labs(color = "Species") +
  theme_dendro() +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )
p2


pdf("figures/pop_gen/phylogeny/shared_buscos_ML_tree_polar.pdf", height = 10, width = 12)
p2
dev.off()


# with Nd
leaves = dendro %>% dendextend::get_nodes_attr("leaf")
leaf_height = dendro %>% dendextend::get_leaves_attr("height")


#######################################
#dendro = as.dendrogram.phylo(fit)
ddata <- ggdendro::dendro_data(dendro, type = "rectangle")

ddata$labels$leaf_height = leaf_height
colnames(dend.labels.new)[c(1,8)] = c("label", "label.new")
ddata$labels = full_join(ddata$labels, dend.labels.new)
ddata$labels

p1 = ggplot() +
    geom_segment(data = ddata$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_point(data = ddata$labels, aes(x = x, y = leaf_height, color = spp), size = 1) +
    #scale_color_manual(values = c25) +
    scale_shape_manual(values = c(18, 16)) +
    coord_flip() +
    # coord_polar() +
    scale_y_reverse() +
    labs(color = "Species") +
    theme_dendro() +
    theme(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position = c(0.8,0.75)
    )
p1

pdf("figures/pop_gen/phylogeny/shared_buscos_ML_tree.with_Nd.pdf", height = 12, width = 8)
p1
dev.off()

#############################################
#############################################
#
#Trying ape style
plot(fit)
plot(fit, type = "unrooted", show.tip = F)
annot = left_join(data.frame(label = fit$tip.label), sample_metadata)
annot
annot$color = vector(mode = "character", length = nrow(annot))
annot[annot$spp == "Nf", "color"] = "blue"
annot[annot$spp == "Nc", "color"] = "red"
annot[annot$spp == "Nd", "color"] = "green"

plot(fit, type = "unrooted", show.tip = F)
ape::tiplabels(annot$spp, bg=annot$color,
          cex=.5)

#root on Neco
tre2 = root(fit, out = 113)
tre2 = ladderize(tre2)
plot(tre2, , show.tip = F)
tiplabels(annot$spp, bg=annot$color,
          cex=.5)

#root on Nedi
tre2 = root(fit, out = 12)
tre2 = ladderize(tre2)
plot(tre2, , show.tip = F)
tiplabels(annot$spp, bg=annot$color,
          cex=.5)

#this is starting to look like we can just rotate the tree to orient Nc correctly
# I think we should try playing with figTree and also doing all assembly and BUSCO extraction