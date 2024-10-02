library(vegan)
library(dplyr)
library(ggplot2)
#library(pals)

capscale(sample_metadata.scaled ~ 1, distance = "euclidean")
sample_metadata = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")

#col_pal = pals::alphabet(n = 18)
scale_color_manual(values = c25) +
