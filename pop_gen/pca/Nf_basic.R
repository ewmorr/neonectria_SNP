library(vegan)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(FactoMineR)
source("library/ggplot_theme.txt")

#Nf
dist.Nf = read.table("data/Nf/final_tables/rm_dups/FINAL_snp.mac_ge2.LD.pca_analyses.dist", header = F) 
dist.ID.Nf = read.table("data/Nf/final_tables/rm_dups/FINAL_snp.mac_ge2.LD.pca_analyses.dist.id", header = F)
rownames(dist.Nf) = dist.ID.Nf[,1]
colnames(dist.Nf) = dist.ID.Nf[,2]
dist.Nf = dist.Nf %>% as.dist()
dist.Nf

#if we are going to remove early samples we should do so before performing PCA
Nf.pca = capscale(dist.Nf ~ 1, distance = "euclidean")
str(Nf.pca)
Nf.eig_vals = Nf.pca$CA$eig/sum(Nf.pca$CA$eig)
plot(Nf.eig_vals)
Nf.pca_scores = data.frame(Nf.pca$CA$u[,1:6])
Nf.pca_scores$Sequence_label = rownames(Nf.pca_scores)
Nf.clust = hclust(d = dist.Nf, method = "ward.D")
plot(Nf.clust)

Nf.hcpc = HCPC(Nf.pca_scores[,1:3], nb.clust = -1)
plot(Nf.hcpc)
Nd.hcpc = HCPC(Nd.pca_scores[,1:3], nb.clust = -1)
#essentially the largest change in inertia gain is where we stop. See section 3.2 here http://factominer.free.fr/more/HCPC_husson_josse.pdf
#
Nd.hcpc = HCPC(Nf.pca_scores[,1:6], nb.clust = 0)

Nf.sample_metadata = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
Nf.pca_scores.metadata = left_join(Nf.pca_scores, Nf.sample_metadata %>% select(Sequence_label, Site, duration_infection, state, collection_period))

ggplot(Nf.pca_scores.metadata,  aes(MDS1, MDS2, color = state, shape = collection_period)) +
    geom_point() + 
    scale_color_manual(values = c25) +
    labs(color = "State/region") +
    my_gg_theme.def_size

ggplot(Nf.pca_scores.metadata %>% filter(state != "VA" & state != "NH.CW" & state != "WV"),  aes(MDS1, MDS2, color = state)) +
    geom_point() + 
    scale_color_manual(values = c25) +
    labs(color = "State/region") +
    my_gg_theme.def_size

ggplot(Nf.pca_scores.metadata,  aes(MDS1, MDS3, color = state)) +
    geom_point() + 
    scale_color_manual(values = c25) +
    labs(color = "State/region") +
    my_gg_theme.def_size

ggplot(Nf.pca_scores.metadata,  aes(MDS2, MDS3, color = state)) +
    geom_point() + 
    scale_color_manual(values = c25) +
    labs(color = "State/region") +
    my_gg_theme.def_size

##################
#pretty plots

p1 = ggplot(
    data.frame(rank = 1:length(Nf.eig_vals), val = Nf.eig_vals*100),
    aes(x = rank, y = val)
) +
    geom_point(shape = 1) +
    my_gg_theme.def_size +
    scale_x_continuous(limits = c(1, NA), breaks = c(1,30,60,90)) +
    labs(y = "Eigen value (% variance)", x = "PCA axis", title = "a") +
    theme(
        plot.title = element_text(hjust = -0.1, margin = margin(b = -10))
    )

p2 = ggplot(Nf.pca_scores.metadata, aes(MDS1, MDS2, color = state)) +
    geom_point() + 
    scale_color_manual(values = c25) +
    labs(
        color = "State/region", 
        x = paste0("PCA1 (", round(Nf.eig_vals[1],3)*100, "% variance)"),
        y = paste0("PCA2 (", round(Nf.eig_vals[2],3)*100, "% variance)"),
        title = "c"
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.1, margin = margin(b = -10))
    )
p2

p2.zoom = ggplot(Nf.pca_scores.metadata %>% 
            filter(state != "VA" & state != "NH.CW" & state != "WV"), 
        aes(MDS1, MDS2, color = state)
    ) +
    geom_point() + 
    scale_color_manual(values = c25) +
    labs(
        color = "State/region", 
        x = paste0("PCA1 (", round(Nf.eig_vals[1],3)*100, "% variance)"),
        y = paste0("PCA2 (", round(Nf.eig_vals[2],3)*100, "% variance)")
    ) +
    my_gg_theme.def_size
p2.zoom

p3 = ggplot(Nf.pca_scores.metadata, aes(MDS1, MDS3, color = state)) +
    geom_point() + 
    scale_color_manual(values = c25) +
    labs(
        color = "State/region", 
        x = paste0("PCA1 (", round(Nf.eig_vals[1],3)*100, "% variance)"),
        y = paste0("PCA2 (", round(Nf.eig_vals[2],3)*100, "% variance)"),
        title = "c"
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.1, margin = margin(b = -10))
    )
p3


#Nd
dist.Nd = read.table("data/Nd/final_tables/rm_dups/FINAL_snp.mac_ge2.LD.pca_analyses.dist", header = F) 
dist.ID.Nd = read.table("data/Nd/final_tables/rm_dups/FINAL_snp.mac_ge2.LD.pca_analyses.dist.id", header = F)
rownames(dist.Nd) = dist.ID.Nd[,1]
colnames(dist.Nd) = dist.ID.Nd[,2]
dist.Nd = dist.Nd %>% as.dist()
dist.Nd

Nd.pca = capscale(dist.Nd ~ 1, distance = "euclidean") #tested bray-curtis, it doesn't seem to make any difference in the level of clustering
str(Nd.pca)
Nd.eig_vals = Nd.pca$CA$eig/sum(Nd.pca$CA$eig)
plot(Nd.eig_vals)
Nd.pca_scores = data.frame(Nd.pca$CA$u[,1:6])
Nd.pca_scores$Sequence_label = rownames(Nd.pca_scores)
Nd.clust = hclust(d = dist.Nd, method = "ward.D")
plot(Nd.clust)

Nd.hcpc = HCPC(Nd.pca_scores[,1:3], nb.clust = -1)
plot(Nd.hcpc)

Nd.sample_metadata = read.csv("data/sample_metadata/Nd_filtered.lat_lon_dur_inf.csv")
Nd.pca_scores.metadata = left_join(Nd.pca_scores, Nd.sample_metadata %>% select(Sequence_label, Site, duration_infection, state, collection_period))

Nd.pca_scores.metadata %>% group_by(state, collection_period) %>% summarize(n = n())

ggplot(Nd.pca_scores.metadata,  aes(MDS1, MDS2, color = state, shape = collection_period)) +
    geom_point() + 
    scale_color_manual(values = c25) +
    labs(color = "State/region") +
    my_gg_theme.def_size

ggplot(Nd.pca_scores.metadata,  aes(MDS1, MDS2, color = state)) +
    geom_point() + 
    scale_color_manual(values = c25) +
    labs(color = "State/region") +
    my_gg_theme.def_size

ggplot(Nd.pca_scores.metadata,  aes(MDS1, MDS3, color = state)) +
    geom_point() + 
    scale_color_manual(values = c25) +
    labs(color = "State/region") +
    my_gg_theme.def_size

ggplot(Nd.pca_scores.metadata,  aes(MDS2, MDS3, color = state)) +
    geom_point() + 
    scale_color_manual(values = c25) +
    labs(color = "State/region") +
    my_gg_theme.def_size

##################
#pretty plots

p3 = ggplot(
    data.frame(rank = 1:length(Nd.eig_vals), val = Nd.eig_vals*100),
    aes(x = rank, y = val)
) +
    geom_point(shape = 1) +
    my_gg_theme.def_size +
    scale_x_continuous(limits = c(1, NA), breaks = c(1,10,20,30)) +
    labs(y = "Eigen value (% variance)", x = "PCA axis", title = "b") +
    theme(
        plot.title = element_text(hjust = -0.1, margin = margin(b = -10))
    )
p3

p4 = ggplot(Nd.pca_scores.metadata, aes(MDS1, MDS2, color = state)) +
    geom_point() + 
    scale_color_manual(values = c25) +
    labs(
        color = "State/region", 
        x = paste0("PCA1 (", round(Nd.eig_vals[1],3)*100, "% variance)"),
        y = paste0("PCA2 (", round(Nd.eig_vals[2],3)*100, "% variance)"),
        title = "d"
    ) +
    my_gg_theme.def_size +
    theme(
        plot.title = element_text(hjust = -0.1, margin = margin(b = -10))
    )
p4


#Nc
dist.Nc = read.table("data/Nc/final_tables/FINAL_snp.mac_ge2.LD.pca_analyses.dist", header = F) 
dist.ID.Nc = read.table("data/Nc/final_tables/FINAL_snp.mac_ge2.LD.pca_analyses.dist.id", header = F)
rownames(dist.Nc) = dist.ID.Nc[,1]
colnames(dist.Nc) = dist.ID.Nc[,2]
dist.Nc = dist.Nc %>% as.dist()
dist.Nc

Nc.pca = capscale(dist.Nc ~ 1, distance = "euclidean")
str(Nc.pca)
Nc.eig_vals = Nc.pca$CA$eig/sum(Nc.pca$CA$eig)
plot(Nc.eig_vals)
Nc.pca_scores = data.frame(scores(Nc.pca)$sites)
Nc.pca_scores$Sequence_label = rownames(Nc.pca_scores)

Nc.sample_metadata = read.csv("data/sample_metadata/Nc_canton_loc_date.lat_lon.csv")
Nc.pca_scores.metadata = left_join(Nc.pca_scores, Nc.sample_metadata %>% select(Sequence_label, Canton))

ggplot(Nc.pca_scores.metadata,  aes(MDS1, MDS2, color = Canton)) +
    geom_point() + 
    scale_color_manual(values = c25) +
    labs(color = "Canton") +
    my_gg_theme.def_size

##################
#pretty plots

p5 = ggplot(
    data.frame(rank = 1:length(Nc.eig_vals), val = Nc.eig_vals*100),
    aes(x = rank, y = val)
) +
    geom_point(shape = 1) +
    my_gg_theme.def_size +
    scale_x_continuous(limits = c(1, NA), breaks = c(1,30,60,90)) +
    labs(y = "Eigen value (% variance)", x = "PCA axis")
p5

p6 = ggplot(Nc.pca_scores.metadata, aes(MDS1, MDS2, color = Canton)) +
    geom_point() + 
    scale_color_manual(values = c25) +
    labs(
        color = "Canton", 
        x = paste0("PCA1 (", round(Nc.eig_vals[1],3)*100, "% variance)"),
        y = paste0("PCA2 (", round(Nc.eig_vals[2],3)*100, "% variance)")
    ) +
    my_gg_theme.def_size
p6


############################
#print plots

pdf("figures/pop_gen/pca/Nf-Nd.scree-PCA.pdf", width = 10, height = 7)
grid.arrange(p1,p3,p2,p4,ncol = 2, heights = c(0.4,0.6))
dev.off()
