#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --job-name="fastHaN"
#SBATCH --output=fastHaN.out
#SBATCH --partition=shared
##SBATCH --exclude=node117,node118

module purge


cd ~/Nd_hapNet
~/repo/fastHaN/fastHaN_linux original_tcs -i FINAL_snp.IBD_analyses.min1.phy -t 24 -o tcs_out

cd ~/Nf_hapNet
~/repo/fastHaN/fastHaN_linux original_tcs -i FINAL_snp.IBD_analyses.min1.phy -t 24 -o tcs_out
