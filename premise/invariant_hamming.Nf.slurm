#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="adegenet"
#SBATCH --output=adegenet.Nf_tigs.out
#SBATCH --partition=shared
#SBATCH --mem=512000

module purge
module load anaconda/colsa
conda activate adegenet

cd ~/nf_nd_ham_dist/Nf/tigs

#this script is a huge memory hog with the larger dataset. Ended up running on a big mem node and it pulled up 350G memory (Nf dataset is ~42Mbp acorss 115 seqs/individuals) (and took a while to run)
# if want to keep using this, should put the file read in a while loop at the very least

#perl ~/repo/neonectria_SNP/library/snp_table2fasta.pl FINAL_invariant.IBD_analyses.table.snps_only 5 > FINAL_invariant.snps_only.for_IBD.fasta

#sed -i 's/NGT//' FINAL_invariant.snps_only.for_IBD.fasta
#Rscript ~/repo/neonectria_SNP/pop_gen/IBD/IBD.individuals.calculate_hamming.invariant.Nf.r

# for individual tigs

for i in tig*snps_only
do(
    tigName=${i%.table.snps_only}
    echo $tigName
    
    perl ~/repo/neonectria_SNP/library/snp_table2fasta.pl $i 5 > $tigName.fasta
    sed -i 's/NGT//' $tigName.fasta
    Rscript ~/repo/neonectria_SNP/pop_gen/IBD/IBD.individuals.calculate_hamming.invariant.Nf.r $tigName.fasta $tigName
)done
