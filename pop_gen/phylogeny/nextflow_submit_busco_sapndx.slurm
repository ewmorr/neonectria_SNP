#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="nextflow"
#SBATCH --output=nextflow.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=24
#SBATCH --exclude=node105,node106,node107,node108,node117,node118,node119,node167
##SBATCH --exclude=node\[101-112,117-125\]
##https://premise.sr.unh.edu/slurm.html#srun-excludenode117node118-hostname
##excluding a bunch bc they are often used by priority users and we don't
## need the special resources
## Except the node range syntax isn't working so we will list it out

cd ~/Nf_Nd_Nc_buscos_SPANDx/

module purge
module load anaconda/colsa

conda activate spandx
source ~/.bashrc

nextflow run ~/SPANDx_git_clone/
