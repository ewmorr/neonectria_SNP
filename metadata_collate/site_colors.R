library(dplyr)
library(ggplot2)
source("library/ggplot_theme.txt")

metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
metadata.Nd = read.csv("data/sample_metadata/Nd_filtered.lat_lon_dur_inf.csv")

all_metadata = rbind(metadata.Nf, metadata.Nd)
nrow(all_metadata)
all_metadata %>% pull(state) %>% unique() %>% length()

Nf.states = metadata.Nf %>% select(state) %>% unique()
nrow(Nf.states)


Nd.states = metadata.Nd %>% select(state) %>% unique()
nrow(Nd.states)

shared = intersect(Nf.states$state, Nd.states$state)
Nf_unique = Nf.states$state[!Nf.states$state %in% shared]
Nd_unique = Nd.states$state[!Nd.states$state %in% shared]

Nf.states$state = c(shared, Nf_unique)
Nd.states$state = c(shared, Nd_unique)

Nf.states$colors = c25[1:nrow(Nf.states)]
Nd.states$colors = c25[1:nrow(Nd.states)]

nrow(Nf.states)
nrow(Nd.states)

write.csv(Nf.states, "data/sample_metadata/Nf_site_colors.csv", quote = F, row.names = F)
write.csv(Nd.states, "data/sample_metadata/Nd_site_colors.csv", quote = F, row.names = F)


