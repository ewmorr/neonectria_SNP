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
#this is the same method used for choosing the number of clusters...
q_vals = vector(length = length(Nf.eig_vals)-2, mode = "numeric")
for(i in 2:length(Nf.eig_vals)-1){
    q_vals[i-1] = (Nf.eig_vals[i]+Nf.eig_vals[i+1])/(Nf.eig_vals[i-1]-Nf.eig_vals[i])
}
plot(q_vals)
q_vals
#keeping in mind that the index is -1 to the axis, 3 axes is the right number
#i.e., the minimum Q

Nf.pca_scores = data.frame(Nf.pca$CA$u[,1:6])
Nf.pca_scores$Sequence_label = rownames(Nf.pca_scores)
#perform the HCPC
Nf.hcpc = HCPC(Nf.pca_scores[,1:3], nb.clust = 0)
plot(Nf.hcpc)
Nf.hcpc$desc.ind
Nf.hcpc$data.clust
Nf.pca_scores = Nf.hcpc$data.clust
Nf.pca_scores$Sequence_label = rownames(Nf.pca_scores)
Nf.pca_scores.metadata = left_join(Nf.pca_scores, Nf.sample_metadata %>% select(Sequence_label, Site, duration_infection, state, collection_period))

plot(Nf.hcpc, choice = "bar")
plot(Nf.hcpc, choice = "tree")

p1 = ggplot(
    data.frame(rank = 1:length(Nf.eig_vals), val = Nf.eig_vals*100),
    aes(x = rank, y = val)
) +
    geom_point(shape = 1) +
    my_gg_theme.def_size +
    scale_x_continuous(limits = c(1, NA), breaks = c(1,30,60,90)) +
    labs(y = "Eigen value (% variance)", x = "PCA axis", title = "a") +
    theme(
        plot.title = element_text(hjust = -0.12, vjust = -2)
    )
p1
p2 = ggplot(Nf.pca_scores.metadata, aes(MDS1, MDS2, color = state)) +
    geom_point() + 
    scale_color_manual(values = c25) +
    labs(
        color = "State/region", 
        x = paste0("PCA1 (", round(Nf.eig_vals[1],3)*100, "% variance)"),
        y = paste0("PCA2 (", round(Nf.eig_vals[2],3)*100, "% variance)"),
        title = "b"
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.15, vjust = -2)
    )
p2


p3 = ggplot(Nf.pca_scores.metadata, aes(MDS1, MDS2, color = clust)) +
    geom_point() + 
    scale_color_manual(values = c25) +
    geom_mark_ellipse(expand = 0, aes(fill=clust), show.legend = F) +
    labs(
        color = "Cluster", 
        x = paste0("PCA1 (", round(Nf.eig_vals[1],3)*100, "% variance)"),
        y = paste0("PCA2 (", round(Nf.eig_vals[2],3)*100, "% variance)"),
        title = "c"
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.15, vjust = -2)
    )
p3

pdf("figures/pop_gen/pca/Nf.clusters.pdf", width = 16, height = 4)
grid.arrange(p1,p2,p3, ncol = 3, widths = c(0.25,0.415,0.335))
dev.off()

par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(Nf.hcpc, choice = "tree", rect = F)
pdf("figures/pop_gen/pca/Nf.clustering.pdf", width = par('din')[1], par('din')[2])
plot(Nf.hcpc, choice = "tree", rect = F)
dev.off()
par(mar = c(5.1, 4.1, 4.1, 2.1))

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
#this is the same method used for choosing the number of clusters...
q_vals = vector(length = length(Nd.eig_vals)-2, mode = "numeric")
for(i in 2:length(Nd.eig_vals)-1){
    q_vals[i-1] = (Nd.eig_vals[i]+Nd.eig_vals[i+1])/(Nd.eig_vals[i-1]-Nd.eig_vals[i])
}
plot(q_vals)
q_vals
#keeping in mind that the index is -1 to the axis, 2 axes is the right number
#i.e., the minimum Q

Nd.pca_scores = data.frame(Nd.pca$CA$u[,1:6])
Nd.pca_scores$Sequence_label = rownames(Nd.pca_scores)
#perform the HCPC
Nd.hcpc = HCPC(Nd.pca_scores[,1:3], nb.clust = 0)
plot(Nd.hcpc)
Nd.hcpc$desc.ind
Nd.hcpc$data.clust
Nd.pca_scores = Nd.hcpc$data.clust
Nd.pca_scores$Sequence_label = rownames(Nd.pca_scores)
Nd.pca_scores.metadata = left_join(Nd.pca_scores, Nd.sample_metadata %>% select(Sequence_label, Site, duration_infection, state, collection_period))

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
p2 = ggplot(Nd.pca_scores.metadata, aes(MDS1, MDS2, color = state)) +
    geom_point() + 
    scale_color_manual(values = c25) +
    labs(
        color = "State/region", 
        x = paste0("PCA1 (", round(Nd.eig_vals[1],3)*100, "% variance)"),
        y = paste0("PCA2 (", round(Nd.eig_vals[2],3)*100, "% variance)"),
        title = "b"
    ) +
    guides(color=guide_legend(ncol=2)) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.15, vjust = -2)
    )
p2


p3 = ggplot(Nd.pca_scores.metadata, aes(MDS1, MDS2, color = clust)) +
    geom_point() + 
    scale_color_manual(values = c25) +
    geom_mark_ellipse(expand = 0, aes(fill=clust), show.legend = F) +
    labs(
        color = "Cluster", 
        x = paste0("PCA1 (", round(Nd.eig_vals[1],3)*100, "% variance)"),
        y = paste0("PCA2 (", round(Nd.eig_vals[2],3)*100, "% variance)"),
        title = "c"
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.15, vjust = -2)
    )
p3

pdf("figures/pop_gen/pca/Nd.clusters.pdf", width = 16, height = 4)
grid.arrange(p1,p2,p3, ncol = 3, widths = c(0.25,0.415,0.335))
dev.off()

par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(Nd.hcpc, choice = "tree", rect = F)
pdf("figures/pop_gen/pca/Nd.clustering.pdf", width = par('din')[1], par('din')[2])
plot(Nd.hcpc, choice = "tree", rect = F)
dev.off()
par(mar = c(5.1, 4.1, 4.1, 2.1))
