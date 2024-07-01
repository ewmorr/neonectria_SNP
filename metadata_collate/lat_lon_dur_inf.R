library(dplyr)

sample_data = read.csv("data/sample_metadata/sample_ID_mapping_all_samples_05092024.csv")
head(sample_data)
site_data = read.csv("data/sample_metadata/site_info.csv")
head(site_data)

all_metadata = left_join(sample_data, site_data, by = "Site")
nrow(all_metadata) == nrow(sample_data)

colnames(all_metadata)
write.csv(
    all_metadata, 
    "data/sample_metadata/lat_lon_dur_inf.csv",
    quote = F,
    row.names = F
)
