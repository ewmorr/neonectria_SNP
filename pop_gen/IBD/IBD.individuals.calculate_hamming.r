library(adegenet)

#########
#########
#We were originally using adegenet to calculate populations level distances
# We would like to instead compute pairwise distances between inds to deal
# with low sample size of Nd
# The genind format can not handle poly alleles (only bialleles)
# The most basic method to calculate distance is using dist() euclidean
# but this fails for polyalleles because different character states will be 
# treated as greater than simple REF:ALT (e.g., 1-0 versus 2-0 vs 2-1).
# Can calculate pairwise Hamming distance (character difs) using 
# phangorn::dist.dna with pairwise deletion and using the nt fasta written
# from the VCF (as used for the phylogeny). Use model = "raw" to compute 
# ratio/proportion of differeing sites and then multiply by the number of
# sites for count based Hamming (which is suitable for poppr::poppr.amova, e.g.)
# 
# ape::dist.dna(x = fasta, model = "raw", pairwise.deletion = T)
# OR
# phangorn::dist.hamming(x = fasta, ratio = T, exclude = "pairwise")
#
#We will use the adegenet/ape version as they give the same answer and the
# AMOVA test we will use is poppr based on adegenet

# read SNP multi seq aligment fasta
#seqs.Nf = phangorn::read.phyDat(
#    file = "data/Nf/final_tables/rm_dups/FINAL_snp.snps_only.for_phylogeny.fasta",
#    format = "fasta",
#    type = "DNA"
#)

#dist.Nf = phangorn::dist.hamming(x = seqs.Nf, ratio = T, exclude = "pairwise")
#head(dist.Nf * 997132)
#[1] 78968.56 78801.49 74581.56 73785.58 75296.76 74680.12

# read SNP multi seq aligment fasta
dnaBin.Nf = adegenet::fasta2DNAbin(file = "data/Nf/final_tables/rm_dups/FINAL_snp.snps_only.for_phylogeny.fasta")
#997,088 nucleotides 

dist.raw.Nf = ape::dist.dna(x = dnaBin.Nf, model = "raw", pairwise.deletion = T)
#head(dist.raw.Nf * 997088)
# 78966.65 78799.58 74579.76 73783.83 75294.95 74678.31

saveRDS(dist.raw.Nf, "data/Nf/IBD/hamming_dist.rds")
rm(dnaBin.Nf)
gc()

#Nd
dnaBin.Nd = adegenet::fasta2DNAbin(file = "data/Nd/final_tables/rm_dups/FINAL_snp.snps_only.for_phylogeny.fasta")
# 1,413,296 nt

dist.raw.Nd = ape::dist.dna(x = dnaBin.Nd, model = "raw", pairwise.deletion = T)

saveRDS(dist.raw.Nd, "data/Nd/IBD/hamming_dist.rds")
rm(dnaBin.Nd)
gc()

#Nc
dnaBin.Nc = adegenet::fasta2DNAbin(file = "data/Nc/final_tables/FINAL_snp.snps_only.for_phylogeny.fasta")
# 378,780

dist.raw.Nc = ape::dist.dna(x = dnaBin.Nc, model = "raw", pairwise.deletion = T)

saveRDS(dist.raw.Nc, "data/Nc/IBD/hamming_dist.rds")
rm(dnaBin.Nc)
gc()


##########################
##########################
#Then run invariants

dnaBin.Nf = adegenet::fasta2DNAbin(file = "data/Nf/final_tables/rm_dups/FINAL_invariant.snps_only.for_IBD.fasta")
#x nucleotides 

dist.raw.Nf = ape::dist.dna(x = dnaBin.Nf, model = "raw", pairwise.deletion = T)

saveRDS(dist.raw.Nf, "data/Nf/IBD/hamming_dist.invariant.rds")
rm(dnaBin.Nf)
gc()

#Nd
dnaBin.Nd = adegenet::fasta2DNAbin(file = "data/Nd/final_tables/rm_dups/FINAL_invariant.snps_only.for_IBD.fasta", chunkSize = 100)
#38,317,560 nucleotides 

dist.raw.Nd = ape::dist.dna(x = dnaBin.Nd, model = "raw", pairwise.deletion = T)

saveRDS(dist.raw.Nd, "data/Nd/IBD/hamming_dist.invariant.rds")
rm(dnaBin.Nd)
gc()
