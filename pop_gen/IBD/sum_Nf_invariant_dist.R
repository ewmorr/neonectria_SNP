######################################
#THIS ROUTINE DID NOT WORK THE NF ESTIMATES ARE WAY TOO HIGH
#WE WILL GET RID OF THIS BUT SAVING FOR NOW


#########################################
#Nf, len norm distances and sum over tigs
test_dist = readRDS("data/Nf/IBD/test_tig00000001_pilon_raw.rds")
range(test_dist)

dist_path = "data/Nf/IBD/invariant"
tig_lens = read.table("data/Nf/IBD/invariant/adegenet.tigs_lens.txt", header = F)
tig_names = read.table("data/Nf/IBD/invariant/adegenet.tigs_names.txt", header = F)
dist_files = list.files(path = dist_path, pattern = "raw.rds", full.names = F)
#everytyhing is in the same sort order so we won't sort further

dist_list = list()
for( i in 1:length(dist_files)){
    temp.dist = readRDS(file.path(dist_path, dist_files[i]))
    print(tig_names[i,1])
    print(tig_lens[i,1])
    #print(min(temp.dist) )
    print(min(temp.dist) * tig_lens[i,1])
    dist_list[[tig_names[i,1]]] = temp.dist * tig_lens[i,1]
}
sum(as.matrix(dist_list[[1]]) < 10^6)

sum_dist = dist_list[[1]]
for( i in 2:length(dist_list)){
    sum_dist = sum_dist + dist_list[[i]]
}
head(sum_dist)
sum_dist_kb = sum_dist/1000
saveRDS(sum_dist, file.path(dist_path, "sum_tigs_len_corrected.rds"))
range(sum_dist[!is.na(sum_dist)])



#########################################
#Nd, len norm distance 
dist = readRDS("data/Nd/IBD/invariant/hamming_dist.invariant.rds")
range(dist)
len = read.table("data/Nd/IBD/invariant/len.txt")
len
dist_len = dist*(len$V1) #per kilobase
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
