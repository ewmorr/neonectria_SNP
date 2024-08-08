library(adegenet)
library(ape)
#library(Biostrings)

#Nf.fasta = readBStringSet(file = "~/nf_nd_ham_dist/Nf/FINAL_invariant.snps_only.for_IBD.fasta", format = "fasta", nrec = 10, use.names=TRUE, seqtype="DNA")
#dnaBin = as.DNAbin(Nf.fasta)
#dnaBin = adegenet::fasta2DNAbin(file = "~/nf_nd_ham_dist/Nf/FINAL_invariant.snps_only.for_IBD.fasta", chunkSize = 30) #use high mem node with 100 chunk
#x nucleotides

args = commandArgs(trailingOnly=TRUE)

fastaFile = args[1]
tigName = args[2]
dnaBin = adegenet::fasta2DNAbin(file = fastaFile, chunkSize = 30) #use high mem node with 100 chunk

dist.raw = ape::dist.dna(x = dnaBin, model = "raw", pairwise.deletion = T)
#dist.raw.len = dist.raw*length(matrix(dnaBin[1,]))

saveRDS(dist.raw, paste(tigName, "raw.rds", sep = ""))
#saveRDS(dist.raw, paste(tigName, "len_corrected.rds", sep = ""))
