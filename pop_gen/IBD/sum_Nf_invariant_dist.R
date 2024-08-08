
#########################################
#Nf, len norm distances and sum over tigs
dist_path = "data/Nf/IBD/invariant"
tig_lens = read.table("data/Nf/IBD/invariant/adegenet.tigs_lens.txt", header = F)
tig_names = read.table("data/Nf/IBD/invariant/adegenet.tigs_names.txt", header = F)
dist_files = list.files(path = dist_path, pattern = "raw.rds", full.names = F)
#everytyhing is in the same sort order so we won't sort further

dist_list = list()
for( i in 1:length(dist_files)){
    temp.dist = readRDS(file.path(dist_path, dist_files[i]))
    dist_list[[tig_names[i,1]]] = temp.dist * tig_lens[i,1]
}

sum_dist = dist_list[[1]]
for( i in 2:length(dist_list)){
    sum_dist = sum_dist + dist_list[[i]]
}
head(sum_dist)
sum_dist_kb = sum_dist/1000
saveRDS(sum_dist, file.path(dist_path, "sum_tigs_len_corrected.rds"))
range(sum_dist[!is.na(sum_dist)])

range(dist_list[[2]])

#########################################
#Nd, len norm distance 
dist = readRDS("data/Nd/IBD/invariant/hamming_dist.invariant.rds")
range(dist)
len = read.table("data/Nd/IBD/invariant/len.txt")
len
dist_len = dist*(len$V1/1000) #per kilobase
saveRDS(dist_len, "data/Nd/IBD/invariant/hamming_dist.len_corrected.rds")
range(dist_len)


#########################################
#Nc, len norm distance 
dist = readRDS("data/Nc/IBD/invariant/hamming_dist.invariant.rds")
range(dist)
len = read.table("data/Nc/IBD/invariant/len.txt")
len
dist_len = dist*(len$V1/1000)
saveRDS(dist_len, "data/Nc/IBD/invariant/hamming_dist.len_corrected.rds")
range(dist_len)
