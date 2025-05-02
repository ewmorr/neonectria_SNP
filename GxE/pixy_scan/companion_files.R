library(dplyr)
    
##################################
pca_clust = read.csv("data/sample_metadata/Nf_filtered.HCPC_clust.csv")
sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
pca_clust$HCPC.cluster = paste0("clust_", pca_clust$HCPC.cluster)

pca_clust.2 = pca_clust %>% filter(HCPC.cluster == "clust_2")

write.table(pca_clust.2, "data/sample_metadata/Nf_pixy_pop.cluster2_single_pop.tsv", sep = "\t", col.names = F, row.names = F, quote = F)

left_join(
    pca_clust.2,
    sample_metadata.Nf %>% select(Sequence_label, state)
) %>% filter(
    state != "NC"
) %>% select(
    -state
) -> pca_clust.2.noNC
nrow(pca_clust.2.noNC)
nrow(pca_clust.2)

write.table(pca_clust.2.noNC, "data/sample_metadata/Nf_pixy_pop.cluster2_noNC_single_pop.tsv", sep = "\t", col.names = F, row.names = F, quote = F)

left_join(
    pca_clust.2,
    sample_metadata.Nf %>% select(Sequence_label, state)
) %>% select(
    -HCPC.cluster
) -> pca_clust.2.statePop
nrow(pca_clust.2.statePop)
head(pca_clust.2.statePop)

write.table(pca_clust.2.statePop, "data/sample_metadata/Nf_pixy_pop.cluster2_site_pops.tsv", sep = "\t", col.names = F, row.names = F, quote = F)

left_join(
    pca_clust.2,
    sample_metadata.Nf %>% select(Sequence_label, state)
) %>% filter(
    state != "NC"
) %>% select(
    -HCPC.cluster
) -> pca_clust.2.noNC.statePop
nrow(pca_clust.2.noNC.statePop)
head(pca_clust.2.noNC.statePop)

write.table(pca_clust.2.noNC.statePop, "data/sample_metadata/Nf_pixy_pop.cluster2_noNC_site_pops.tsv", sep = "\t", col.names = F, row.names = F, quote = F)
