library(dplyr)
library(ggplot2)
library(phylogram)
library(ggdendro)
library(dendextend)
#library(zoo)
library(ape)
library(TreeTools)

#############
#Read tree
source("library/ggplot_theme.txt")



tre = treeio::read.tree("data/Fugr1_ref/map_against_Fugr_no_core_extract/core.Fusgr1-neonectria.snps_aln.ph")
tre = treeio::read.tree("data/Fugr1_ref/map_against_Fugr_no_core_extract/fugr_root.nwk")
plot(tre)

tre = RootTree(tre, "Fusgr1_core_snp_pos")
tre = ape::root.phylo(tre, "Fusgr1_core_snp_pos", resolve.root = T)

dendro = as.dendrogram.phylo(tre)
plot(dendro)

#Read metadata
#
sample_metadata = read.csv("data/sample_metadata/core_fugr.lat_lon_dur_inf.csv")
colnames(sample_metadata)
colnames(sample_metadata)[1] = "label" #for join below
#
#########

# first plot

dend.labels.new = left_join(data.frame(label = dendro %>% labels), sample_metadata)
dend.labels.new$label.new = paste0(dend.labels.new$spp, 1:nrow(dend.labels.new))
dend.labels.new$color = vector(mode = "character", length = nrow(dend.labels.new))
dend.labels.new[dend.labels.new$spp == "Nf", "color"] = "#D55E00"
dend.labels.new[dend.labels.new$spp == "Nc", "color"] = "#009E73"
dend.labels.new[dend.labels.new$spp == "Nd", "color"] = "#0072B2"
dend.labels.new[dend.labels.new$spp == "Fg", "color"] = "black"

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
    dendextend::ladderize(., right = T) %>%
    plot(horiz = T)

dendro %>% 
    set("labels", dend.labels.new$label.new) %>%
    set("labels_colors", dend.labels.new$color) %>%
    #prune(Nd.labels) %>%
    #rank_branches %>%
    dendextend::ladderize(., right = T) %>%
    plot(horiz = T)

dendro %>% 
    set("labels", dend.labels.new$label.new) %>%
    set("labels_colors", dend.labels.new$color) %>%
    prune(Nd.labels) %>%
    #rank_branches %>%
    dendextend::ladderize(., right = T) %>%
    circlize_dendrogram() %>%
    plot(horiz = T)

dendro %>% 
    set("labels", dend.labels.new$label.new) %>%
    set("labels_colors", dend.labels.new$color) %>%
    prune(Nd.labels) %>%
    rank_branches %>%
    dendextend::ladderize(., right = T) %>%
    circlize_dendrogram()

###################
#this one is decent
dendro %>% 
    set("labels", dend.labels.new$label.new) %>%
    set("labels_colors", dend.labels.new$color) %>%
    #prune(Nd.labels) %>%
    rank_branches %>%
    dendextend::ladderize(., right = T) %>%
    circlize_dendrogram()
##################
##################

dendro %>% 
    set("labels", dend.labels.new$label.new) %>%
    set("labels_colors", dend.labels.new$color) %>%
#    prune(Nd.labels) %>%
    rank_branches %>%
    dendextend::ladderize(., right = T) %>%
    hang.dendrogram() %>% # this is not good it cuts off the branch len
    circlize_dendrogram()

#####################
#####################
# similar to the one we like above but with points
# instead of text
pdf("figures/pop_gen/phylogeny/fugr_outgroup.polar.ladderize.points.pdf")
dendro %>% 
    set("labels", dend.labels.new$label.new) %>%
    set("leaves_pch", 19) %>%
    set("leaves_col", dend.labels.new$color) %>%
    #    prune(Nd.labels) %>%
    #rank_branches %>%
    #hang.dendrogram() %>%
    dendextend::ladderize(., right = T) %>%
    circlize_dendrogram(labels = F)
dev.off()
#####################
#####################

dendro %>% 
    set("labels", dend.labels.new$label.new) %>%
    set("labels_colors", dend.labels.new$color) %>%
    #set("leaves_pch", 19) %>%
    #set("leaves_col", dend.labels.new$color) %>%
    #prune(Nd.labels) %>%
    rank_branches %>%
    dendextend::ladderize(., right = T) %>%
#    hang.dendrogram() %>%
    plot(horiz = T)
# nope



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
  scale_color_manual(values = cbPalette) +
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

pdf("figures/pop_gen/phylogeny/fugr_outgroup.no_Nd.rectangle.pdf", height = 12, width = 8)
p1
dev.off()


p2 = ggplot() +
  geom_segment(data = ddata$segments, aes(x = x, y = y/100, xend = xend, yend = yend/100)) +
  geom_point(data = ddata$labels, aes(x = x, y = leaf_height/100, color = spp), size = 1) +
  #scale_color_manual(values = c25) +
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


#pdf("figures/pop_gen/phylogeny/shared_buscos_ML_tree_polar.pdf", height = 10, width = 12)
#p2
#dev.off()


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
    scale_color_manual(values = cbPalette) +
    scale_shape_manual(values = c(18, 16)) +
    coord_flip() +
    #coord_polar() +
    scale_y_reverse() +
    labs(color = "Species") +
    theme_dendro() +
    theme(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)#,
        #legend.position = c(0.8,0.75)
    )
p1

pdf("figures/pop_gen/phylogeny/fugr_outgroup.with_Nd.rectangle.pdf", height = 12, width = 8)
p1
dev.off()

p2 = ggplot() +
    geom_segment(data = ddata$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_point(data = ddata$labels, aes(x = x, y = leaf_height, color = spp), size = 1) +
    scale_color_manual(values = cbPalette) +
    scale_shape_manual(values = c(18, 16)) +
    coord_flip() +
    coord_polar() +
    scale_y_reverse() +
    labs(color = "Species") +
    theme_dendro() +
    theme(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)#,
        #legend.position = c(0.8,0.75)
    )
p2

pdf("figures/pop_gen/phylogeny/fugr_outgroup.polar.pdf", height = 12, width = 8)
p2
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