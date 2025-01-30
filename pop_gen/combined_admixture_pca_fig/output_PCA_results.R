library(vegan)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(FactoMineR)
#library(ggforce)
#source("library/ggplot_theme.txt")


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

###########
Nf.pca_scores = Nf.hcpc$data.clust
Nf.pca_scores$Sequence_label = rownames(Nf.pca_scores)
Nf.pca_scores.metadata = left_join(Nf.pca_scores, Nf.sample_metadata %>% select(Sequence_label, Site, duration_infection, state, collection_period, Tree_species))
write.csv(Nf.pca_scores.metadata, "data/pca/Nf_pca.metadata.csv", row.names = F, quote = F)


######################
# Nd
######################
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

###########
Nd.pca_scores = Nd.hcpc$data.clust
Nd.pca_scores$Sequence_label = rownames(Nd.pca_scores)
Nd.pca_scores.metadata = left_join(Nd.pca_scores, Nd.sample_metadata %>% select(Sequence_label, Site, duration_infection, state, collection_period, Tree_species))
write.csv(Nd.pca_scores.metadata, "data/pca/Nd_pca.metadata.csv", row.names = F, quote = F)

