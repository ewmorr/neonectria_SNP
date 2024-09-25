#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="bbduk"
#SBATCH --output=bbduk.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa
conda activate qiime2-2019.4

inDir=/mnt/home/garnas/ewj4/core_genome_fugr_SPANDx
outDir=/mnt/home/garnas/ewj4/neonectria_all_reads_for_analysis_trimmed

adapterFile=/mnt/lustre/software/anaconda/colsa/envs/qiime2-2019.4/opt/bbmap-38.57-0/resources/adapters.fa

for i in $inDir/*R1*.fastq.gz
do

    file1=${i##*/}
    sampleName=${file1%R1*}
    postName=${file1##*R1}
    fwd=$file1
    rev=${sampleName}R2${postName}

srun bbduk.sh in=$inDir/$fwd in2=$inDir/$rev \
out=$outDir/$fwd out2=$outDir/$rev \
ref=$adapterFile ktrim=r k=23 mink=11 hdist=1 tpe tbo \
qtrim=r trimq=10 -minlength=36

done
