library(pegas)
library(vcfR) #we'll read and convert with the vcfR functions

vcf <- read.vcfR("data/Nf_SPANDx_all_seqs/out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.vcf.gz", verbose = FALSE)
ref = read.FASTA("data/ref.fasta")

dna = vcfR2DNAbin(x = vcf, ref.seq = ref, start.pos = 1)
#this memory dumps. May need to precess by chromosome
dna = vcfR2DNAbin(x = vcf)

seg.sites(dna) %>% length
length(dna) #