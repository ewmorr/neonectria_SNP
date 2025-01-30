library(vegan)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(FactoMineR)
library(ggforce)
source("library/ggplot_theme.txt")


#Nd
dist.Nd = read.table("data/Nd/final_tables/rm_dups/FINAL_snp.mac_ge2.LD.pca_analyses.dist", header = F) 
dist.ID.Nd = read.table("data/Nd/final_tables/rm_dups/FINAL_snp.mac_ge2.LD.pca_analyses.dist.id", header = F)
rownames(dist.Nd) = dist.ID.Nd[,1]
colnames(dist.Nd) = dist.ID.Nd[,2]
class(dist.Nd)
Nd.sample_metadata = read.csv("data/sample_metadata/Nd_filtered.lat_lon_dur_inf.csv")
#if we are going to remove early samples we should do so before performing PCA
#rm_ids = Nd.sample_metadata %>% filter(collection_period == "early") %>% pull(Sequence_label)
#dist.Nd[rownames(dist.Nd) %in% rm_ids,] = NA
#dist.Nd[,colnames(dist.Nd) %in% rm_ids] = NA
#sum(is.na(dist.Nd))
#convert to dist
dist.Nd = dist.Nd %>% as.dist()
dist.Nd



#PCa with capscale (could use FactoMineR::PCA but we know this way already)
Nd.pca = capscale(dist.Nd ~ 1, distance = "euclidean", na.action = na.omit)
str(Nd.pca)
Nd.eig_vals = Nd.pca$CA$eig/sum(Nd.pca$CA$eig)
plot(Nd.eig_vals)

Nd.pca_scores = data.frame(Nd.pca$CA$u[,1:6])
Nd.pca_scores$Sequence_label = rownames(Nd.pca_scores)
#perform the HCPC
Nd.hcpc = HCPC(Nd.pca_scores[,1:3], nb.clust = 0)
plot(Nd.hcpc)
Nd.hcpc$desc.ind
Nd.hcpc$data.clust
#write for metadata
write.csv(
    data.frame(Sequence_label = row.names(Nd.hcpc$data.clust), HCPC.cluster = Nd.hcpc$data.clust$clust),
    "data/sample_metadata/Nd_filtered.HCPC_clust.csv",
    row.names = F,
    quote = F
)
###########
Nd.pca_scores = Nd.hcpc$data.clust
Nd.pca_scores$Sequence_label = rownames(Nd.pca_scores)
Nd.pca_scores.metadata = left_join(Nd.pca_scores, Nd.sample_metadata %>% select(Sequence_label, Site, duration_infection, state, collection_period, Tree_species))

tree_adonis = adonis2(dist.Nd ~ Tree_species, data = Nd.pca_scores.metadata)
summary(tree_adonis)
tree_adonis
# Tree_species  8 1.7239e+10 0.28248 1.0334  0.073 .
adonis2(dist.Nd ~ state, data = Nd.pca_scores.metadata)
# state    17 3.6689e+10 0.60119 1.0641  0.019 *
adonis2(dist.Nd ~ state*Tree_species, data = Nd.pca_scores.metadata)
#state        17 3.6689e+10 0.60119 1.0673  0.017 *
#Tree_species  1 2.0962e+09 0.03435 1.0367  0.098 .
early_collections = Nd.pca_scores.metadata %>% filter(collection_period == "early") %>% pull(Sequence_label)
dist.Nd.no_early = as.matrix(dist.Nd)
dist.Nd.no_early[rownames(dist.Nd) == early_collections,] = NA
dist.Nd.no_early[,colnames(dist.Nd) == early_collections] = NA

adonis2(dist.Nd.no_early ~ duration_infection, data = Nd.pca_scores.metadata)
#duration_infection  1 2.8149e+09 0.04612 1.3539  0.001 ***
adonis2(dist.Nd.no_early ~ Tree_species, data = Nd.pca_scores.metadata)
#Tree_species  8 1.7239e+10 0.28248 1.0334  0.055 .
adonis2(dist.Nd.no_early ~ Tree_species+duration_infection, data = Nd.pca_scores.metadata)
#
adonis2(dist.Nd.no_early ~ Tree_species+duration_infection+state, data = Nd.pca_scores.metadata)
#Run test of all with partial (by = "margin")
anova(dbrda(dist.Nd.no_early ~ state+Tree_species+duration_infection, data=Nd.pca_scores.metadata, distance = "bray"), by="margin")
#state               9 666671848 1.0624  0.024 *
#    Tree_species        1  72281175 1.0367  0.114  
#duration_infection  0         0   -Inf         
anova(dbrda(dist.Nd.no_early ~ state+duration_infection+Tree_species, data=Nd.pca_scores.metadata, distance = "bray"), by="margin")
anova(dbrda(dist.Nd.no_early ~ duration_infection+Tree_species, data=Nd.pca_scores.metadata, distance = "bray"), by="margin")
#duration_infection  1   76285988 1.0642  0.107
#Tree_species        8  573673530 1.0004  0.398
anova(dbrda(dist.Nd.no_early ~ state+Tree_species, data=Nd.pca_scores.metadata, distance = "bray"), by="margin")
anova(dbrda(dist.Nd.no_early ~ duration_infection+state, data=Nd.pca_scores.metadata, distance = "bray"), by="margin")
anova(dbrda(dist.Nd.no_early ~ duration_infection, data=Nd.pca_scores.metadata, distance = "bray"), by="terms")
# duration_infection  1   97064012 1.3539  0.001 ***
anova(dbrda(dist.Nd.no_early ~ state+duration_infection+Tree_species, data=Nd.pca_scores.metadata, distance = "bray"), by="terms")
# state is always significant (margins) and explains 60% of variance alone (terms)


plot(Nd.hcpc, choice = "bar")
plot(Nd.hcpc, choice = "tree")

p1 = ggplot(
    data.frame(rank = 1:length(Nd.eig_vals), val = Nd.eig_vals*100),
    aes(x = rank, y = val)
) +
    geom_point(shape = 1) +
    my_gg_theme.def_size +
    scale_x_continuous(limits = c(1, NA), breaks = c(1,10,20,30)) +
    labs(y = "Eigen value (% variance)", x = "PCA axis", title = "a") +
    theme(
        plot.title = element_text(hjust = -0.12, vjust = -2)
    )
p1
p2 = ggplot(Nd.pca_scores.metadata, 
        aes(MDS1, MDS2, 
            color = state, 
            shape = factor(collection_period, 
                levels = c("early", "modern"), 
                labels = c("1960s", "contemporary")
            )
        )
    ) +
    geom_point(size = 2, position = position_jitter(width = 0.015)) + 
    scale_color_manual(values = c25) +
    scale_shape_manual(values = c(17,16)) +
    labs(
        color = "State/region", 
        shape = "Collection period",
        x = paste0("PCA1 (", round(Nd.eig_vals[1],3)*100, "% variance)"),
        y = paste0("PCA2 (", round(Nd.eig_vals[2],3)*100, "% variance)"),
        title = "a"
    ) +
    guides(color=guide_legend(ncol=2)) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.19, vjust = -2) #hjust -0.2 for main, -0.17 for supp
    ) 
p2

p3 = ggplot(Nd.pca_scores.metadata, 
        aes(MDS1, MDS2, 
            color = Tree_species, 
            shape = factor(collection_period, 
                levels = c("early", "modern"), 
                labels = c("1960s", "contemporary")
            )
        )
    ) +
    geom_point(size = 2, position = position_jitter(width = 0.015)) + 
    scale_color_brewer(palette = "Set1") +
    scale_shape_manual(values = c(17,16)) +
    guides(shape = "none") + #comment for supp plot
    labs(
        shape = "Collection period",
        color = "Tree species", 
        x = paste0("PCA1 (", round(Nd.eig_vals[1],3)*100, "% variance)"),
        y = paste0("PCA2 (", round(Nd.eig_vals[2],3)*100, "% variance)"),
        title = "b" #b for main, c for supp
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.19, vjust = -2) #hjust -0.19 for main, -0.17 for supp
    )
p3


p4 = ggplot(Nd.pca_scores.metadata, aes(MDS1, MDS2, color = clust)) +
    geom_point() + 
    scale_color_manual(values = c25) +
    geom_mark_ellipse(expand = 0, aes(fill=clust), show.legend = F) +
    labs(
        color = "Cluster", 
        x = paste0("PCA1 (", round(Nd.eig_vals[1],3)*100, "% variance)"),
        y = paste0("PCA2 (", round(Nd.eig_vals[2],3)*100, "% variance)"),
        title = "c" #c for main, d for supp 
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.19, vjust = -2) #hjust -0.19 for main, -0.17 for supp
    )
p4

p5 = ggplot(Nd.pca_scores.metadata %>% filter(collection_period == "modern"), 
            aes(MDS1, MDS2, 
                color = duration_infection
            )
) +
    geom_point(size = 2, position = position_jitter(width = 0.015)) + 
    labs(
        color = "Infestation\nduration (yrs)", 
        x = paste0("PCA1 (", round(Nd.eig_vals[1],3)*100, "% variance)"),
        y = paste0("PCA2 (", round(Nd.eig_vals[2],3)*100, "% variance)"),
        title = "b"
    ) +
    scale_color_gradient2(high = "#2166ac", low = "#b2182b", mid = "#d1e5f0", midpoint = 50) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.17, vjust = -2)
    )
p5

pdf("figures/pop_gen/pca/Nd.tree_sp.pdf", width = 16, height = 4)
grid.arrange(p2,p3,p4, ncol = 3, widths = c(0.37, 0.325, 0.305))
dev.off()

library(gtable)
gp2 = ggplotGrob(p2)
gp3 = ggplotGrob(p3)
gp5 = ggplotGrob(p5)
gp4 = ggplotGrob(p4)

gp23 = rbind(gp2,gp3)
gp54 = rbind(gp5,gp4)
gpall = cbind(gp23, gp54)

pdf("figures/pop_gen/pca/Nd.tree_sp.inf_dur.pdf", width = 12, height = 8)
#grid.arrange(p2,p5,p3,p4, ncol = 2, widths = c(0.55, 0.45))
plot(gpall)
dev.off()

