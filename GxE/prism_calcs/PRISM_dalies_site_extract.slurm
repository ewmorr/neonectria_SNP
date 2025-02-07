#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="PRISM"
#SBATCH --output=PRISM.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa

conda activate PRISM

srun Rscript ~/repo/neonectria_SNP/GxE/PRISM_calcs/tmin_tmax_ppt_dailys_site_extract.r
