#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --job-name="fastHaN"
#SBATCH --output=fastHaN.out
#SBATCH --partition=shared
#SBATCH --exclude=node117,node118

module purge


cd ~/Nd_hapNet
~/repo/fastHaN/fastHaN_linux mjn -i ~/Nd_hapNet/FINAL_snp.IBD_analyses.min1.phy -t 24 -o ~/Nd_hapNet/tcs_out

cd ~/Nf_hapNet
~/repo/fastHaN/fastHaN_linux mjn -i ~/Nf_hapNet/FINAL_snp.IBD_analyses.min1.phy -t 24 -o ~/Nf_hapNet/tcs_out
