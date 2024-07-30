library(adegenet)


dnaBin = adegenet::fasta2DNAbin(file = "~/nf_nd_ham_dist/Nf/FINAL_invariant.snps_only.for_IBD.fasta", chunkSize = 25)
#x nucleotides

dist.raw = ape::dist.dna(x = dnaBin, model = "raw", pairwise.deletion = T)

saveRDS(dist.raw, "~/nf_nd_ham_dist/Nf/hamming_dist.invariant.rds")
