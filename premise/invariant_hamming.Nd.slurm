#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="adegenet"
#SBATCH --output=adegenet.Nd.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa
conda activate adegenet

cd ~/nf_nd_ham_dist/Nd/
#perl ~/repo/neonectria_SNP/library/snp_table2fasta.pl FINAL_invariant.IBD_analyses.table.snps_only 5 > FINAL_invariant.snps_only.for_IBD.fasta
#sed -i 's/NGT//' FINAL_invariant.snps_only.for_IBD.fasta
Rscript ~/repo/neonectria_SNP/pop_gen/IBD/IBD.individuals.calculate_hamming.invariant.Nd.r

