source("R_scripts/make_site_metadata.r")
source("R_scripts/make_structure_indfile.r")

low_n_states = ( (ind_file %>% group_by(state.name) %>% summarize(n = n()) ) %>% filter(n < 4) )$state.name

retain_samps = sample_metadata %>% 
                    filter(!state.name %in% low_n_states & sample %in% ind_file$sample)  %>%
    select(sample)

write.table(retain_samps, file = "data/Nf_SPANDx_all_seqs/retain_samples.txt", quote = F, col.names = F, row.names = F)
