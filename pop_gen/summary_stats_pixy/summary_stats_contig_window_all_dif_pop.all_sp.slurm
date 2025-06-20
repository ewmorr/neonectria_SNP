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

outdir=contig_all_dif_pop

indir=$HOME/Nf_pixy
pixy --stats dxy \
    --vcf $indir/FINAL_invariant.nuclear.vcf.gz \
    --populations $indir/Nf.all_dif_pixy_pop.tsv \
    --bed_file contigs.bed \
    --output_folder $indir/$outdir \
    --n_cores 24

indir=$HOME/Nd_pixy
pixy --stats dxy \
    --vcf $indir/FINAL_invariant.nuclear.vcf.gz \
    --populations $indir/Nd.all_dif_pixy_pop.tsv \
    --bed_file contigs.bed \
    --output_folder $indir/$outdir \
    --n_cores 24

indir=$HOME/Nc_pixy
pixy --stats dxy \
    --vcf $indir/FINAL_invariant.nuclear.vcf.gz \
    --populations $indir/Nc.all_dif_pixy_pop.tsv \
    --bed_file contigs.bed \
    --output_folder $indir/$outdir \
    --n_cores 24

