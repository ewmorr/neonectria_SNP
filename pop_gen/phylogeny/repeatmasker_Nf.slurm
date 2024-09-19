#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="repeatmasker"
#SBATCH --cpus-per-task=12
#SBATCH --output=repeatmasker.Nf.out
#SBATCH --partition=shared
#SBATCH --exclude=node117,node118

module purge
module load linuxbrew/colsa


cd Nf_SPANDx_all_seqs

RepeatMasker -s -pa 12 -lib ~/RepeatMaskerLib/fngrep.ref ref.fasta
