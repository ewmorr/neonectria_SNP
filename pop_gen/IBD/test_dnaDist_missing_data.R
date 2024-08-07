library(ape)

x <- structure(c("As", "N", "hyphen", "G", 
                 "AAAAAAAAAA", "AAAAAAAAAN", "AAAAAAAAA-", "AAAAAAAAAG"), .Dim = c(4L, 2L))
x
y <- t(sapply(strsplit(x[,2],""), as.character))
y
rownames(y) <- x[,1]
yBin = as.DNAbin(y)

ape::dist.dna(yBin, model = "raw", pairwise.deletion = T)


tig_lens = read.table("data/Nf/final_tables/rm_dups/invariant_table_tigs/lens.txt")
sum(tig_lens$V1)

length(str(yBin))
printlen(yBin)

length(matrix(yBin[1,]))

tig_lens = read.table("data/Nf/final_tables/rm_dups/invariant_table_tigs/lens_adegenet.txt")
sum(tig_lens)

tig_lens = read.table("data/Nf/final_tables/rm_dups/invariant_table_tigs_retry/lens.txt")
sum(tig_lens$V1)
