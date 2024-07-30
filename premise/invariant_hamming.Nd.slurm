#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="adegenet"
#SBATCH --output=adegenet.Nd.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa
conda activate adegenet

Rscript ~/repo/neonectria_SNP/pop_gen/IBD/IBD.individuals.calculate_hamming.invariant.Nd.r

