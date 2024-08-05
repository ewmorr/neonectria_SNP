library(adegenet)
#library(Biostrings)

#Nf.fasta = readBStringSet(file = "~/nf_nd_ham_dist/Nf/FINAL_invariant.snps_only.for_IBD.fasta", format = "fasta", nrec = 10, use.names=TRUE, seqtype="DNA")
#dnaBin = as.DNAbin(Nf.fasta)

dnaBin = adegenet::fasta2DNAbin(file = "~/nf_nd_ham_dist/Nf/FINAL_invariant.snps_only.for_IBD.fasta", chunkSize = 30) #use high mem node with 100 chunk
#x nucleotides

dist.raw = ape::dist.dna(x = dnaBin, model = "raw", pairwise.deletion = T)

saveRDS(dist.raw, "~/nf_nd_ham_dist/Nf/hamming_dist.invariant.fast_mem.rds")
