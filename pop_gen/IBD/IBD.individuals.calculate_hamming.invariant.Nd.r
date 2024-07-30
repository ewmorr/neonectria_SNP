library(adegenet)


dnaBin = adegenet::fasta2DNAbin(file = "~/nf_nd_ham_dist/Nd/FINAL_invariant.snps_only.for_IBD.fasta")
#x nucleotides

dist.raw = ape::dist.dna(x = dnaBin, model = "raw", pairwise.deletion = T)

saveRDS(dist.raw, "~/nf_nd_ham_dist/Nd/hamming_dist.invariant.rds")
