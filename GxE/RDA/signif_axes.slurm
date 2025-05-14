#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="cca.full"
#SBATCH --output=anova.cca.full.out
#SBATCH --partition=shared
#SBATCH --exclude=node117,node118
#SBATCH --cpus-per-task=24

module purge
module load anaconda/colsa
conda activate r-4.2

# first arg is the input file and the second is the output file which will be written to the same dir
Rscript ~/repo/neonectria_SNP/GxE/RDA/signif_axes.R $HOME/GxE_RDA/no_NC.RDA_env_full.rds signif_axis.env_full.rds

