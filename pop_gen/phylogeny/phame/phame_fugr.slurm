#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="phame"
#SBATCH --output=phame.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --exclude=node105,node106,node107,node108,node117,node118,node119,node167
##SBATCH --exclude=node\[101-112,117-125\]
##https://premise.sr.unh.edu/slurm.html#srun-excludenode117node118-hostname
##excluding a bunch bc they are often used by priority users and we don't
## need the special resources
## Except the node range syntax isn't working so we will list it out

cd ~/fugr_neonectria_phame/

module purge
module load anaconda/colsa

conda activate phame_new

phame phame.ctl

