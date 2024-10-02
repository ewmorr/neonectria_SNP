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

pdf("figures/pop_gen/admixture/Nf.admixture.site_facets.pdf", width = 10, height = 5)
p1.mod
dev.off()
