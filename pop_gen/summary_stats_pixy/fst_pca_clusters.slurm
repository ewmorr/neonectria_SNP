#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name="pixy"
#SBATCH --output=pixy.out
#SBATCH --partition=shared
#SBATCH --exclude=node105,node106,node107,node108,node117,node118,node119,node167

module purge
module load anaconda/colsa
conda activate pixy

cd ~/Nf_pixy/

pixy --stats fst \
    --vcf FINAL_invariant.nuclear.vcf.gz \
    --populations Nf_pixy_pops.HCPC_clust.txt \
    --bed_file contigs.bed \
    --n_cores 12

mv pixy_*.txt whole_contig_clusters_pop

pixy --stats pi fst dxy \
    --vcf FINAL_invariant.nuclear.vcf.gz \
    --populations Nf_pixy_pops.HCPC_clust.txt \
    --window_size 10000 \
    --n_cores 12

mv pixy_*.txt window_10Kb_clusters_pop/
