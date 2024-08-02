library(adegenet)


dnaBin = adegenet::fasta2DNAbin(file = "~/nf_nd_ham_dist/Nf/FINAL_invariant.snps_only.for_IBD.fasta", chunkSize = 125) #use high mem node with 100 chunk
#x nucleotides

dist.raw = ape::dist.dna(x = dnaBin, model = "raw", pairwise.deletion = T)

saveRDS(dist.raw, "~/nf_nd_ham_dist/Nf/hamming_dist.invariant.fast_mem.rds")
