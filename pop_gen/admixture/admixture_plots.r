library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
source("library/ggplot_theme.txt")

#Nf metadata and ancestry as list
Nf.fam_info = read.table("data/Nf/final_tables/rm_dups/FINAL_snp.admixture.nosex", header = F)
colnames(Nf.fam_info)[1] = "Sequence_label"
Nf.sample_metadata = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
Nf.sample_metadata = left_join(data.frame(Sequence_label = Nf.fam_info[,1]), Nf.sample_metadata)

Nf.K_list = list()
for(i in 2:14){
    file_name = paste0("data/Nf/admixture/FINAL_snp.admixture.",i,".Q")
    Nf.K_list[[i]] = read.table(file_name, header = F)
    colnames(Nf.K_list[[i]]) = paste0("Q", 1:ncol(Nf.K_list[[i]]))
    Nf.K_list[[i]]$Sequence_label = Nf.fam_info$Sequence_label
    Nf.K_list[[i]]$K = i
    Nf.K_list[[i]] = Nf.K_list[[i]]  %>% 
        pivot_longer(cols = starts_with("Q"), names_to = "ancestor", values_to = "Q")
}
Nf.K_df = bind_rows(Nf.K_list)
Nf.K_df = left_join(
    Nf.K_df, 
    Nf.sample_metadata %>% 
        select(Sequence_label, state, duration_infection, collection_period)
)
#Nd metadata and ancestry as list
Nd.fam_info = read.table("data/Nd/final_tables/rm_dups/FINAL_snp.admixture.nosex", header = F)
colnames(Nd.fam_info)[1] = "Sequence_label"
Nd.sample_metadata = read.csv("data/sample_metadata/Nd_filtered.lat_lon_dur_inf.csv")
Nd.sample_metadata = left_join(data.frame(Sequence_label = Nd.fam_info[,1]), Nd.sample_metadata)

Nd.K_list = list()
for(i in 2:14){
    file_name = paste0("data/Nd/admixture/FINAL_snp.admixture.",i,".Q")
    Nd.K_list[[i]] = read.table(file_name, header = F)
    colnames(Nd.K_list[[i]]) = paste0("Q", 1:ncol(Nd.K_list[[i]]))
    Nd.K_list[[i]]$Sequence_label = Nd.fam_info$Sequence_label
    Nd.K_list[[i]]$K = i
    Nd.K_list[[i]] = Nd.K_list[[i]] %>% 
        pivot_longer(cols = starts_with("Q"), names_to = "ancestor", values_to = "Q")
}
Nd.K_df = bind_rows(Nd.K_list)
Nd.K_df = left_join(
    Nd.K_df, 
    Nd.sample_metadata %>% 
        select(Sequence_label, state, duration_infection, collection_period)
)




########################################
# plots
# 

#sort by infection duration then by state
Nf.try_order = Nf.sample_metadata[with(Nf.sample_metadata, order(duration_infection, state)),"Sequence_label"] 
Nd.try_order = Nd.sample_metadata[with(Nd.sample_metadata, order(duration_infection, state)),"Sequence_label"] 

K_number = 9
p1 = ggplot(Nf.K_df %>% filter(collection_period == "modern" & K <= K_number), 
       aes(
           x = factor(Sequence_label, levels = Nf.try_order), 
           y = Q, 
           fill = ancestor
        )
    ) +
    geom_bar(stat = "identity") +
    facet_wrap(
        ~factor(
            K, 
            levels = 2:K_number, 
            labels = paste0("K = ", 2:K_number) 
        ),
        nrow = K_number,
        switch = "y"
    ) +
    scale_fill_brewer(palette = "Set1", guide = "none") +
    my_gg_theme.def_size +
    labs(
        x = "Individuals (ordered by infestation duration and site)",
        y = "Ancestry"
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside"
    )
p1

p2 = ggplot(Nd.K_df %>% filter(collection_period == "modern" & K <= K_number), 
            aes(
                x = factor(Sequence_label, levels = Nd.try_order), 
                y = Q, 
                fill = ancestor
            )
) +
    geom_bar(stat = "identity") +
    facet_wrap(
        ~factor(
            K, 
            levels = 2:K_number, 
            labels = paste0("K = ", 2:K_number) 
        ),
        nrow = K_number,
        switch = "y"
    ) +
    scale_fill_brewer(palette = "Set1", guide = "none") +
    my_gg_theme.def_size +
    labs(
        x = "Individuals (ordered by infestation duration and site)",
        y = "Ancestry"
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside"
    )
p2

pdf("figures/pop_gen/admixture/admixture.pdf", width = 18, height = 7)
grid.arrange(p1,p2, widths = c(0.8, 0.2))
dev.off()


#######################
#k4 with state grouping 
#
Nf.K_df$facet_by_state = paste0(Nf.K_df$state, "--",Nf.K_df$duration_infection, " yrs")
Nf.K_df[Nf.K_df$state %in% c("NB.LUD", "NB.RB"), "facet_by_state"] = "NB 1963 collection" # we lump the NB early collections

Nf.K_df %>% select(state, duration_infection) %>% unique() %>% print(n = Inf)
Nf.K_df %>% select(state, duration_infection, collection_period) %>% unique() %>% print(n = Inf)

Nf.unique_collections = Nf.K_df %>% select(duration_infection, collection_period, facet_by_state) %>% unique() 

state_order = Nf.unique_collections[with(Nf.unique_collections, order(collection_period, duration_infection)) , ] %>% pull(facet_by_state) 

p1.mod = ggplot(Nf.K_df %>% filter(K == 4), 
aes(
    x = Sequence_label, 
    y = Q, 
    fill = ancestor
)
) +
    geom_bar(stat = "identity") +
    facet_wrap(~factor(facet_by_state, levels = rev(state_order)), scales = "free_x", ncol = 7) + 
    scale_fill_brewer(palette = "Set1", guide = "none") +
    my_gg_theme.def_size +
    labs(
        x = "Individuals (unique sites in different facets ordered by infestation duration)",
        y = "Ancestry (K = 4)"
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside"
    )
p1.mod

#Nd
Nd.K_df$facet_by_state = paste0(Nd.K_df$state, "--",Nd.K_df$duration_infection, " yrs")
Nd.K_df %>% group_by(facet_by_state, collection_period, state) %>% summarize(n = n())
#NB.NE -- 1964
#QC.C -- 1959
#QC.E -- 1965

Nd.K_df[Nd.K_df$state %in% c("NB.NE", "QC.C", "QC.E") & Nd.K_df$collection_period == "early", "facet_by_state"] = "NB/QC 1959-1965" # we lump the early collections

Nd.K_df %>% select(state, duration_infection) %>% unique() %>% print(n = Inf)
Nd.K_df %>% select(state, duration_infection, collection_period) %>% unique() %>% print(n = Inf)

Nd.unique_collections = Nd.K_df %>% select(duration_infection, collection_period, facet_by_state) %>% unique() 

state_order = Nd.unique_collections[with(Nd.unique_collections, order(collection_period, duration_infection)) , ] %>% pull(facet_by_state) %>% unique

p2.mod = ggplot(Nd.K_df %>% filter(K == 5), 
                aes(
                    x = Sequence_label, 
                    y = Q, 
                    fill = ancestor
                )
) +
    geom_bar(stat = "identity") +
    facet_wrap(~factor(facet_by_state, levels = rev(state_order)), scales = "free_x", ncol = 8) + 
    scale_fill_brewer(palette = "Set1", guide = "none") +
    my_gg_theme.def_size +
    labs(
        x = "Individuals (unique sites in different facets ordered by infestation duration)",
        y = "Ancestry (K = 5)"
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside"
    )
p2.mod

    

pdf("figures/pop_gen/admixture/Nf.admixture.site_facets.pdf", width = 10, height = 5)
p1.mod
dev.off()

pdf("figures/pop_gen/admixture/Nd.admixture.site_facets.pdf", width = 10, height = 5)
p2.mod
dev.off()



###############################
# K 2-5; sort by inf dur then
# by state as above

Nf.K_df$facet_by_state = paste0(Nf.K_df$state, "--",Nf.K_df$duration_infection, " yrs")
Nf.K_df[Nf.K_df$state %in% c("NB.LUD", "NB.RB"), "facet_by_state"] = "NB 1963 collection" # we lump the NB early collections

Nf.K_df %>% select(state, duration_infection) %>% unique() %>% print(n = Inf)
Nf.K_df %>% select(state, duration_infection, collection_period) %>% unique() %>% print(n = Inf)

Nf.unique_collections = Nf.K_df %>% select(duration_infection, collection_period, facet_by_state) %>% unique() 

state_order = Nf.unique_collections[with(Nf.unique_collections, order(collection_period, duration_infection)) , ] %>% pull(facet_by_state) 

Nf_clust = read.csv("data/sample_metadata/Nf_filtered.HCPC_clust.csv")
Nf.K_df.clust = left_join(Nf.K_df, Nf_clust)

#sort individuals within facets by cluster
#the facets handle the state/dur inf sorting
p1.mod = ggplot(Nf.K_df.clust %>% filter(K == 4), 
                aes(
                    x = reorder(Sequence_label, HCPC.cluster), 
                    #x = Sequence_label,
                    y = Q, 
                    fill = ancestor
                )
) +
    geom_bar(stat = "identity") +
    ggforce::facet_row(
        ~factor(facet_by_state, levels = rev(state_order)), 
        scales = "free_x", 
        space = "free",
        strip.position = "left"
    ) + 
    scale_fill_brewer(palette = "Set1", guide = "none") +
    my_gg_theme.def_size +
    labs(
        x = "Individuals (unique sites in different facets ordered by infestation duration)",
        y = "Ancestry (K = 4)"
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside"
    )
p1.mod


###############################
# K 2-5; sort by  state
# facet by HCPC cluster

Nf.K_df$facet_by_state = paste0(Nf.K_df$state, "--",Nf.K_df$duration_infection, " yrs")
Nf.K_df[Nf.K_df$state %in% c("NB.LUD", "NB.RB"), "facet_by_state"] = "NB 1963 collection" # we lump the NB early collections

Nf.K_df %>% select(state, duration_infection) %>% unique() %>% print(n = Inf)
Nf.K_df %>% select(state, duration_infection, collection_period) %>% unique() %>% print(n = Inf)

Nf.unique_collections = Nf.K_df %>% select(duration_infection, collection_period, facet_by_state) %>% unique() 

state_order = Nf.unique_collections[with(Nf.unique_collections, order(collection_period, duration_infection)) , ] %>% pull(facet_by_state) 

# read and join cluster data
Nf_clust = read.csv("data/sample_metadata/Nf_filtered.HCPC_clust.csv")
Nf.K_df.clust = left_join(Nf.K_df, Nf_clust)



#sort individuals within facets by cluster
#the facets handle the state/dur inf sorting
p1.mod = ggplot(Nf.K_df.clust %>% filter(K == 4), 
                aes(
                    x = reorder(Sequence_label, HCPC.cluster), 
                    #x = Sequence_label,
                    y = Q, 
                    fill = ancestor
                )
) +
    geom_bar(stat = "identity") +
    ggforce::facet_row(
        ~factor(facet_by_state, levels = rev(state_order)), 
        scales = "free_x", 
        space = "free",
        strip.position = "left"
    ) + 
    scale_fill_brewer(palette = "Set1", guide = "none") +
    my_gg_theme.def_size +
    labs(
        x = "Individuals (unique sites in different facets ordered by infestation duration)",
        y = "Ancestry (K = 4)"
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside"
    )
p1.mod

# grid on HCPC clust and K
Nf.try_order = Nf.sample_metadata[with(Nf.sample_metadata, order(duration_infection, state)),"Sequence_label"] 

#remaking clust as factor

p1.mod = ggplot(Nf.K_df.clust %>% filter(K >= 2 & K <= 5), 
                aes(
                    x = factor(Sequence_label, levels = Nf.try_order), 
                    y = Q, 
                    fill = ancestor
                )
) +
    geom_bar(stat = "identity") +
    facet_grid(
        factor(K,
               levels = c(2,3,4,5),
               labels = paste("K = ", 2:5)
        )~factor(HCPC.cluster, 
                 levels = c(1,2,3,4), # 4,1,3,2
                 labels = c("1: WV", "2: mixed sites including VA and NH.CW individuals", "3: NH.CW", "4: VA")
        ), 
        scales = "free_x", 
        space = "free_x",
        switch = "both"
    ) + 
    scale_fill_brewer(palette = "Set1", guide = "none") +
    my_gg_theme.def_size +
    labs(
        x = "HCPC cluster (bars represent individuals)",
        y = "Ancestry",
        title = "a"
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = -0.02, vjust = -2)
    )
p1.mod

#nd
#
Nd.K_df$facet_by_state = paste0(Nd.K_df$state, "--",Nd.K_df$duration_infection, " yrs")
Nd.K_df %>% group_by(facet_by_state, collection_period, state) %>% summarize(n = n())
#NB.NE -- 1964
#QC.C -- 1959
#QC.E -- 1965

Nd.K_df[Nd.K_df$state %in% c("NB.NE", "QC.C", "QC.E") & Nd.K_df$collection_period == "early", "facet_by_state"] = "NB/QC 1959-1965" # we lump the early collections

Nd.K_df %>% select(state, duration_infection) %>% unique() %>% print(n = Inf)
Nd.K_df %>% select(state, duration_infection, collection_period) %>% unique() %>% print(n = Inf)

Nd.unique_collections = Nd.K_df %>% select(duration_infection, collection_period, facet_by_state) %>% unique() 

state_order = Nd.unique_collections[with(Nd.unique_collections, order(collection_period, duration_infection)) , ] %>% pull(facet_by_state) %>% unique

# read and join cluster data

Nd_clust = read.csv("data/sample_metadata/Nd_filtered.HCPC_clust.csv")
Nd.K_df.clust = left_join(Nd.K_df, Nd_clust)

Nd.try_order = Nd.sample_metadata[with(Nd.sample_metadata, order(duration_infection, state)),"Sequence_label"] 


p2.mod = ggplot(Nd.K_df.clust %>% filter(K >= 2 & K <= 5), 
                aes(
                    x = factor(Sequence_label, levels = Nd.try_order), 
                    y = Q, 
                    fill = ancestor
                )
) +
    geom_bar(stat = "identity") +
    facet_grid(
        factor(K,
               levels = c(2,3,4,5),
               labels = paste("K = ", 2:5)
        )~HCPC.cluster, 
        scales = "free_x", 
        space = "free_x",
        switch = "both"
    ) + 
    scale_fill_brewer(palette = "Set1", guide = "none") +
    my_gg_theme.def_size +
    labs(
        x = "HCPC cluster (bars represent individuals)",
        y = "Ancestry",
        title = "b"
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = -0.075, vjust = -2)
    )
p2.mod



pdf("figures/pop_gen/admixture/K2-5.HCPC_clust_facets.pdf", width = 18, height = 7)
grid.arrange(p1.mod,p2.mod, widths = c(0.77, 0.23))
dev.off()
