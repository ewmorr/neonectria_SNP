#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="adegenet"
#SBATCH --output=adegenet.four_samples.out
#SBATCH --partition=shared
#SBATCH --exclude=node105,node106,node107,node108,node117,node118,node119,node167

module purge
module load anaconda/colsa
conda activate R-pop_gen

cd ~/Nf_pop_IBD_11182024/

Rscript ~/repo/neonectria_SNP/pop_gen/IBD/IBD.adegenet_metrics.multiple_subsamples.min_3.PREMISE.r

