#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --job-name="thin_str_th"
#SBATCH --output=structure_threader.out
#SBATCH --partition=shared
#SBATCH --exclude=node105,node106,node107,node108,node117,node118,node119,node167

module purge
module load anaconda/colsa

export PATH=~/.local/bin:$PATH

cd ~/Nf_SPANDx_all_seqs_structure_th_thin_locData/

#main params

structure_threader run -K 9 -R 5 \
-i ~/Nf_SPANDx_all_seqs_structure_th_thin_locData/FINAL_snp.thinned_no_NA.LocData.structure \
-o ~/Nf_SPANDx_all_seqs_structure_th_thin_locData/ \
--params ~/Nf_SPANDx_all_seqs_structure_th_thin_locData/mainparams \
-t 24 -st /mnt/home/garnas/ewj4/.local/bin/structure

