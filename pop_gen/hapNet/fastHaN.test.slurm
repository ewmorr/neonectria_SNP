#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --job-name="fastHaN"
#SBATCH --output=fastHaN.test.out
#SBATCH --partition=shared
#SBATCH --mem=512G
#SBATCH --exclude=node117,node118

module purge

~/repo/fastHaN/fastHaN_linux original_tcs -i ~/repo/fastHaN/TestData/Test1000.phy.gz -t 8 -o Test1000
