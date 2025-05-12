#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="anova.cca"
#SBATCH --output=anova.cca.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=24

module purge
module load anaconda/colsa
conda activate r-4.2

Rscript ~/repo/neonectria_SNP/GxE/RDA/signif_axes.R $HOME/GxE_RDA/no_NC.RDA.rds

