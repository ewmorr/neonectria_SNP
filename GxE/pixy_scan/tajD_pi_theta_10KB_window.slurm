#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name="pixy"
#SBATCH --output=pixy_10kb.out
#SBATCH --partition=shared
#SBATCH --exclude=node105,node106,node107,node108,node117,node118,node119,node167

module purge
module load anaconda/colsa
conda activate pixy2

cd ~/Nf_pixy_scan_gxe/

mkdir cluster2_single_pop

pixy --stats pi watterson_theta tajima_d \
    --vcf FINAL_invariant.nuclear.cluster2.vcf.gz \
    --populations Nf_pixy_pop.cluster2_single_pop.tsv \
    --window_size 10000 \
    --n_cores 12

mv pixy_*.txt cluster2_single_pop/

mkdir cluster2_noNC_single_pop

pixy --stats pi watterson_theta tajima_d \
    --vcf FINAL_invariant.nuclear.cluster2_no_NC.vcf.gz \
    --populations Nf_pixy_pop.cluster2_noNC_single_pop.tsv \
    --window_size 10000 \
    --n_cores 12

mv pixy_*.txt cluster2_noNC_single_pop/

mkdir cluster2_site_pops

pixy --stats pi watterson_theta tajima_d \
    --vcf FINAL_invariant.nuclear.cluster2.vcf.gz \
    --populations Nf_pixy_pop.cluster2_site_pops.tsv \
    --window_size 10000 \
    --n_cores 12

mv pixy_*.txt cluster2_site_pops/

mkdir cluster2_noNC_site_pops

pixy --stats pi watterson_theta tajima_d \
    --vcf FINAL_invariant.nuclear.cluster2.vcf.gz \
    --populations Nf_pixy_pop.cluster2_noNC_site_pops.tsv \
    --window_size 10000 \
    --n_cores 12

mv pixy_*.txt cluster2_noNC_site_pops/
