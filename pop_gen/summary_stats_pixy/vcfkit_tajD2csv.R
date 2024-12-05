library(dplyr)

sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
state_n = sample_metadata.Nf %>% group_by(state) %>% summarise(n = n())
state_n %>% print(n = Inf)
low_n = state_n %>% filter(n < 2)
#tried up to 5
#in our original dataset where we had sig negative relationship
#there were only seven sites
#see the massmyco presentation from 2021

states = sample_metadata.Nf %>% filter(!state %in% low_n$state) %>% pull(state) %>% unique()


taj_df = data.frame(
    state = states,
    TajimaD = vector(mode = "numeric", length = length(states) )
)

for(i in 1:length(states)){
    file_name = paste0("data/Nf/pixy/", states[i], ".tajima.tsv")
    temp = read.table(file_name, header = T)
    taj_df[i,"state"] = states[i]
    taj_df[i,"TajimaD"] = mean(temp$TajimaD)
}
taj_df
sample_metadata.states = sample_metadata.Nf %>% select(state, duration_infection) %>% distinct()
taj.meta = left_join(taj_df, sample_metadata.states)

plot(TajimaD ~ duration_infection, data = taj.meta)
summary(lm(TajimaD ~ duration_infection, data = taj.meta))
#no stat sig
#
write.csv(taj.meta, "data/Nf/pixy/tajD.vcfkit.csv", row.names = F, quote = F)
