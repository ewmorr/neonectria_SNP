#!/usr/bin/bash
#SBATCH --job-name="gunzip"
#SBATCH --output=gunzip.out
#SBATCH --partition=shared
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --array=0-12%6
#SBATCH --exclude=node117,node118

cd ~/fugr_neonectria_phame/workdir

array=($(find . -name "*.fastq.gz" -printf "%f "))
echo ${array[$SLURM_ARRAY_TASK_ID]}
gunzip ${array[$SLURM_ARRAY_TASK_ID]}
