library(vegan)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(FactoMineR)
library(ggforce)
source("library/ggplot_theme.txt")


#Nf
dist.Nf = read.table("data/Nf/final_tables/rm_dups/FINAL_snp.mac_ge2.LD.pca_analyses.dist", header = F) 
dist.ID.Nf = read.table("data/Nf/final_tables/rm_dups/FINAL_snp.mac_ge2.LD.pca_analyses.dist.id", header = F)
rownames(dist.Nf) = dist.ID.Nf[,1]
colnames(dist.Nf) = dist.ID.Nf[,2]
class(dist.Nf)
Nf.sample_metadata = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
#if we are going to remove early samples we should do so before performing PCA
#rm_ids = Nf.sample_metadata %>% filter(collection_period == "early") %>% pull(Sequence_label)
#dist.Nf[rownames(dist.Nf) %in% rm_ids,] = NA
#dist.Nf[,colnames(dist.Nf) %in% rm_ids] = NA
#sum(is.na(dist.Nf))
#convert to dist
dist.Nf = dist.Nf %>% as.dist()
dist.Nf


#PCa with capscale (could use FactoMineR::PCA but we know this way already)
Nf.pca = capscale(dist.Nf ~ 1, distance = "euclidean", na.action = na.omit)
str(Nf.pca)
Nf.eig_vals = Nf.pca$CA$eig/sum(Nf.pca$CA$eig)
plot(Nf.eig_vals)

Nf.pca_scores = data.frame(Nf.pca$CA$u[,1:6])
Nf.pca_scores$Sequence_label = rownames(Nf.pca_scores)
#perform the HCPC
Nf.hcpc = HCPC(Nf.pca_scores[,1:3], nb.clust = 0)
plot(Nf.hcpc)
Nf.hcpc$desc.ind
Nf.hcpc$data.clust
#write for metadata
write.csv(
    data.frame(Sequence_label = row.names(Nf.hcpc$data.clust), HCPC.cluster = Nf.hcpc$data.clust$clust),
    "data/sample_metadata/Nf_filtered.HCPC_clust.csv",
    row.names = F,
    quote = F
)
###########
Nf.pca_scores = Nf.hcpc$data.clust
Nf.pca_scores$Sequence_label = rownames(Nf.pca_scores)
Nf.pca_scores.metadata = left_join(Nf.pca_scores, Nf.sample_metadata %>% select(Sequence_label, Site, duration_infection, state, collection_period, Tree_species))

tree_adonis = adonis2(dist.Nf ~ duration_infection, data = Nf.pca_scores.metadata)
summary(tree_adonis)
tree_adonis
# duration_infection   1 9.8235e+09 0.02165 2.5004  0.001 ***
adonis2(dist.Nf ~ state, data = Nf.pca_scores.metadata)
# state     22 1.3951e+11 0.30745 1.8565  0.001 ***
adonis2(dist.Nf ~ state*duration_infection, data = Nf.pca_scores.metadata)
#same result
early_collections = Nf.pca_scores.metadata %>% filter(collection_period == "early") %>% pull(Sequence_label)
dist.Nf.no_early = as.matrix(dist.Nf)
dist.Nf.no_early[rownames(dist.Nf) == early_collections,] = NA
dist.Nf.no_early[,colnames(dist.Nf) == early_collections] = NA

adonis2(dist.Nf.no_early ~ duration_infection, data = Nf.pca_scores.metadata)
# duration_infection   1 9.8235e+09 0.02165 2.5004  0.001 ***
#
#Run test of all with partial (by = "margin")
anova(dbrda(dist.Nf.no_early ~ state+duration_infection, data=Nf.pca_scores.metadata, distance = "bray"), by="margin")
#state              21 1137622175 1.808  0.001 ***
#duration_infection  0          0   Inf           
anova(dbrda(dist.Nf.no_early ~ duration_infection+state, data=Nf.pca_scores.metadata, distance = "bray"), by="margin")
#duration_infection  0          0  -Inf           
#state              21 1137622175 1.808  0.001 ***

anova(dbrda(dist.Nf.no_early ~ duration_infection, data=Nf.pca_scores.metadata, distance = "bray"), by="terms")
# duration_infection   1   86171030 2.5004  0.001 ***
anova(dbrda(dist.Nf.no_early ~ duration_infection+state, data=Nf.pca_scores.metadata, distance = "bray"), by="terms")
#duration_infection  1   86171030 2.8759  0.001 ***
#state              21 1137622175 1.8080  0.001 ***
#    Residual           92 2756622681                  

################
################
#Plots

p1 = ggplot(
    data.frame(rank = 1:length(Nf.eig_vals), val = Nf.eig_vals*100),
    aes(x = rank, y = val)
) +
    geom_point(shape = 1) +
    my_gg_theme.def_size +
    #scale_x_continuous(limits = c(1, NA), breaks = c(1,10,20,30)) +
    labs(y = "Eigen value (% variance)", x = "PCA axis", title = "a") +
    theme(
        plot.title = element_text(hjust = -0.12, vjust = -2)
    )
p1
p2 = ggplot(Nf.pca_scores.metadata, 
        aes(MDS1, MDS2, 
            color = state, 
            shape = factor(collection_period, 
                levels = c("early", "modern"), 
                labels = c("1960s", "contemporary")
            )
        )
    ) +
    geom_point(size = 2) + 
    scale_color_manual(values = c25) +
    scale_shape_manual(values = c(17,16)) +
    labs(
        color = "State/region", 
        shape = "Collection period",
        x = paste0("PCA1 (", round(Nf.eig_vals[1],3)*100, "% variance)"),
        y = paste0("PCA2 (", round(Nf.eig_vals[2],3)*100, "% variance)"),
        title = "a"
    ) +
    guides(color=guide_legend(ncol=2)) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.15, vjust = -2) #hjust -0.2 for main, -0.17 for supp
    ) 
p2
p2.mod = ggplot(Nf.pca_scores.metadata %>% filter(!state %in% c("NH.CW", "WV", "VA")), 
            aes(MDS1, MDS2, 
                color = state, 
                shape = factor(collection_period, 
                               levels = c("early", "modern"), 
                               labels = c("1960s", "contemporary")
                )
            )
) +
    geom_point(size = 2) + 
    scale_color_manual(values = c25) +
    scale_shape_manual(values = c(17,16)) +
    labs(
        color = "State/region", 
        shape = "Collection period",
        x = paste0("PCA1 (", round(Nf.eig_vals[1],3)*100, "% variance)"),
        y = paste0("PCA2 (", round(Nf.eig_vals[2],3)*100, "% variance)"),
        title = "a"
    ) +
    guides(color=guide_legend(ncol=2)) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.15, vjust = -2) #hjust -0.2 for main, -0.17 for supp
    ) 
p2.mod

p3 = ggplot(Nf.pca_scores.metadata %>% filter(collection_period == "modern"), 
            aes(MDS1, MDS2, 
                color = duration_infection
            )
) +
    geom_point(size = 2) + 
    labs(
        color = "Infestation\nduration (yrs)", 
        x = paste0("PCA1 (", round(Nf.eig_vals[1],3)*100, "% variance)"),
        y = paste0("PCA2 (", round(Nf.eig_vals[2],3)*100, "% variance)"),
        title = "b"
    ) +
    scale_color_gradient2(high = "#2166ac", low = "#b2182b", mid = "#d1e5f0", midpoint = 45) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.15, vjust = -2)
    )
p3

p3.mod = ggplot(Nf.pca_scores.metadata %>% filter(collection_period == "modern" & !state %in% c("NH.CW", "WV", "VA")), 
                aes(MDS1, MDS2, 
                    color = duration_infection
                )
) +
    geom_point(size = 2) + 
    labs(
        color = "Infestation\nduration (yrs)", 
        x = paste0("PCA1 (", round(Nf.eig_vals[1],3)*100, "% variance)"),
        y = paste0("PCA2 (", round(Nf.eig_vals[2],3)*100, "% variance)"),
        title = "b"
    ) +
    scale_color_gradient2(high = "#2166ac", low = "#b2182b", mid = "#d1e5f0", midpoint = 45) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.15, vjust = -2)
    )
p3.mod

p4 = ggplot(Nf.pca_scores.metadata, aes(MDS1, MDS2, color = clust)) +
    geom_point() + 
    scale_color_manual(values = c25) +
    geom_mark_ellipse(expand = 0, aes(fill=clust), show.legend = F) +
    labs(
        color = "Cluster", 
        x = paste0("PCA1 (", round(Nf.eig_vals[1],3)*100, "% variance)"),
        y = paste0("PCA2 (", round(Nf.eig_vals[2],3)*100, "% variance)"),
        title = "b" #c for main, d for supp, b no scree
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.15, vjust = -2) #hjust -0.19 for main, -0.17 for supp
    )
p4


pdf("figures/pop_gen/pca/Nf.durInf.pdf", width = 19, height = 5)
grid.arrange(p2,p3,p4, ncol = 3, widths = c(0.37, 0.325, 0.305))
dev.off()

pdf("figures/pop_gen/pca/Nf.clusters_no_scree.pdf", width = 12.75, height = 4.75)
grid.arrange(p2,p4, ncol = 2, widths = c(0.55, 0.45))
dev.off()
