library(vcfR)

# NF

#filtered VCF
vcf <- read.vcfR("data/Nf/out.vcf", verbose = FALSE)

#just pulling in dp to get sample IDs

sample_ids <- extract.gt(vcf, element='DP', as.numeric=TRUE) %>% colnames()

#sample_ids
first_set.ids = paste0("NG", 1:99)
second_set.ids = paste0("NG", 101:163)
third_set.ids = c(
    "NG103",
    "NG111",
    "NG112",
    "NG114",
    "NG116",
    "NG121",
    "NG155",
    "NG160",
    "NG117",
    "NG163"
)
second_set.ids = second_set.ids[!second_set.ids %in% third_set.ids]
fourth_set.ids = paste0("NG", 170:196)


lib_one = first_set.ids[first_set.ids %in% sample_ids]
lib_two = second_set.ids[second_set.ids %in% sample_ids]
lib_three = third_set.ids[third_set.ids %in% sample_ids]
lib_four = fourth_set.ids[fourth_set.ids %in% sample_ids]

write.table(lib_one, "data/Nf/lib_one_ids.txt", quote = F, row.names = F, col.names = F)
write.table(lib_two, "data/Nf/lib_two_ids.txt", quote = F, row.names = F, col.names = F)
write.table(lib_three, "data/Nf/lib_three_ids.txt", quote = F, row.names = F, col.names = F)
write.table(lib_four, "data/Nf/lib_four_ids.txt", quote = F, row.names = F, col.names = F)

# ND
#filtered VCF
vcf <- read.vcfR("data/Nd/out.vcf", verbose = FALSE)

#just pulling in dp to get sample IDs

sample_ids <- extract.gt(vcf, element='DP', as.numeric=TRUE) %>% colnames()

#sample_ids
first_set.ids = paste0("NG", 1:99)
second_set.ids = paste0("NG", 101:163)
third_set.ids = c(
    "NG103",
    "NG111",
    "NG112",
    "NG114",
    "NG116",
    "NG121",
    "NG155",
    "NG160",
    "NG117",
    "NG163"
)
second_set.ids = second_set.ids[!second_set.ids %in% third_set.ids]
fourth_set.ids = paste0("NG", 170:196)


lib_one = first_set.ids[first_set.ids %in% sample_ids]
lib_two = second_set.ids[second_set.ids %in% sample_ids]
lib_three = third_set.ids[third_set.ids %in% sample_ids]
lib_four = fourth_set.ids[fourth_set.ids %in% sample_ids]

write.table(lib_one, "data/Nd/lib_one_ids.txt", quote = F, row.names = F, col.names = F)
write.table(lib_two, "data/Nd/lib_two_ids.txt", quote = F, row.names = F, col.names = F)
write.table(lib_three, "data/Nd/lib_three_ids.txt", quote = F, row.names = F, col.names = F)
write.table(lib_four, "data/Nd/lib_four_ids.txt", quote = F, row.names = F, col.names = F)


####################
# Shared BUSCOs

#filtered VCF
vcf <- read.vcfR("data/shared_buscos/out.vcf", verbose = FALSE)

#just pulling in dp to get sample IDs

sample_ids <- extract.gt(vcf, element='DP', as.numeric=TRUE) %>% colnames()

#sample_ids
first_set.ids = paste0("NG", 1:99)
second_set.ids = paste0("NG", 101:163)
third_set.ids = c(
    "NG103",
    "NG111",
    "NG112",
    "NG114",
    "NG116",
    "NG121",
    "NG155",
    "NG160",
    "NG117",
    "NG163"
)
second_set.ids = second_set.ids[!second_set.ids %in% third_set.ids]
fourth_set.ids = paste0("NG", 170:196)
fifth_set.ids = paste0("Neco-", 2:6)


lib_one = first_set.ids[first_set.ids %in% sample_ids]
lib_two = second_set.ids[second_set.ids %in% sample_ids]
lib_three = third_set.ids[third_set.ids %in% sample_ids]
lib_four = fourth_set.ids[fourth_set.ids %in% sample_ids]
lib_five = fifth_set.ids[fifth_set.ids %in% sample_ids]

write.table(lib_one, "data/shared_buscos/lib_one_ids.txt", quote = F, row.names = F, col.names = F)
write.table(lib_two, "data/shared_buscos/lib_two_ids.txt", quote = F, row.names = F, col.names = F)
write.table(lib_three, "data/shared_buscos/lib_three_ids.txt", quote = F, row.names = F, col.names = F)
write.table(lib_four, "data/shared_buscos/lib_four_ids.txt", quote = F, row.names = F, col.names = F)
write.table(lib_five, "data/shared_buscos/lib_five_ids.txt", quote = F, row.names = F, col.names = F)

####################
# Shared BUSCOs

#filtered VCF
vcf <- read.vcfR("data/Fugr1_ref/map_against_Fugr_then_extract_nucmer_core/out.vcf", verbose = FALSE)

#just pulling in dp to get sample IDs

sample_ids <- extract.gt(vcf, element='DP', as.numeric=TRUE) %>% colnames()

#sample_ids
first_set.ids = paste0("NG", 1:99)
second_set.ids = paste0("NG", 101:163)
third_set.ids = c(
    "NG103",
    "NG111",
    "NG112",
    "NG114",
    "NG116",
    "NG121",
    "NG155",
    "NG160",
    "NG117",
    "NG163"
)
second_set.ids = second_set.ids[!second_set.ids %in% third_set.ids]
fourth_set.ids = paste0("NG", 170:196)
fifth_set.ids = paste0("Neco-", 2:6)


lib_one = first_set.ids[first_set.ids %in% sample_ids]
lib_two = second_set.ids[second_set.ids %in% sample_ids]
lib_three = third_set.ids[third_set.ids %in% sample_ids]
lib_four = fourth_set.ids[fourth_set.ids %in% sample_ids]
lib_five = fifth_set.ids[fifth_set.ids %in% sample_ids]

write.table(lib_one, "data/Fugr1_ref/map_against_Fugr_then_extract_nucmer_core/lib_one_ids.txt", quote = F, row.names = F, col.names = F)
write.table(lib_two, "data/Fugr1_ref/map_against_Fugr_then_extract_nucmer_core/lib_two_ids.txt", quote = F, row.names = F, col.names = F)
write.table(lib_three, "data/Fugr1_ref/map_against_Fugr_then_extract_nucmer_core/lib_three_ids.txt", quote = F, row.names = F, col.names = F)
write.table(lib_four, "data/Fugr1_ref/map_against_Fugr_then_extract_nucmer_core/lib_four_ids.txt", quote = F, row.names = F, col.names = F)
write.table(lib_five, "data/Fugr1_ref/map_against_Fugr_then_extract_nucmer_core/lib_five_ids.txt", quote = F, row.names = F, col.names = F)
