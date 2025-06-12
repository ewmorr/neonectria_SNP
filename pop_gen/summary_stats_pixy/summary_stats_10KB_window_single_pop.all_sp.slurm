#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --job-name="pixy"
#SBATCH --output=pixy.out
#SBATCH --partition=shared
#SBATCH --exclude=node105,node106,node107,node108,node117,node118,node119,node167

module purge
module load anaconda/colsa
conda activate pixy2

outdir=window_10kb_single_pop

#indir=$HOME/Nf_pixy
#pixy --stats pi watterson_theta tajima_d \
#    --vcf $indir/FINAL_invariant.nuclear.vcf.gz \
#    --populations $indir/Nf.single_pixy_pop.tsv \
#    --window_size 10000 \
#    --output_folder $indir/$outdir \
#    --n_cores 24

indir=$HOME/Nd_pixy
pixy --stats pi watterson_theta tajima_d \
    --vcf $indir/FINAL_invariant.nuclear.vcf.gz \
    --populations $indir/Nd.single_pixy_pop.tsv \
    --window_size 10000 \
    --output_folder $indir/$outdir \
    --n_cores 24

indir=$HOME/Nc_pixy
pixy --stats pi watterson_theta tajima_d \
    --vcf $indir/FINAL_invariant.nuclear.vcf.gz \
    --populations $indir/Nc.single_pixy_pop.tsv \
    --window_size 10000 \
    --output_folder $indir/$outdir \
    --n_cores 24

