#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="PRISM"
#SBATCH --output=PRISM.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa

conda activate PRISM

srun Rscript ~/repo/neonectria_genome_reseq_10072020/R_scripts/PRISM_dailies_tmin_download.r
