#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --job-name="str_th"
#SBATCH --output=structure_threader.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa

export PATH=~/.local/bin:$PATH

cd ~/Nf_SPANDx_all_seqs_structure_th/

structure_threader run -K 15 -R 5 \
-i ~/Nf_SPANDx_all_seqs/Outputs/Master_vcf/out.filtered.LD_filtered_0.5_10Kb.structure \
-o ~/Nf_SPANDx_all_seqs_structure_th/ \
--params ~/Nf_SPANDx_all_seqs_structure_th/mainparams \
-t 24 -st /mnt/home/garnas/ewj4/.local/bin/structure

