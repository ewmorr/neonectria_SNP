#!/usr/bin/bash
#SBATCH --job-name="gatk"
#SBATCH --output=gatk.out
#SBATCH --partition=shared
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --array=0-146%24
#SBATCH --exclude=node117,node118

module purge
module load anaconda/colsa

conda activate spandx

#mkdir Nd_invariant_sites_GVCF
#mkdir Nd_invariant_sites_GVCF/indv_GVCFs

#change array number appropriately
#find $ref_dir/Outputs/bams/ -name *.dedup.bam
#find $ref_dir/Outputs/bams/ -name *.dedup.bam | wc -l

ref_dir=/mnt/home/garnas/ewj4/core_genome_fugr_SPANDx
out_dir=/mnt/home/garnas/ewj4/core_fugr_invariant_sites_GVCF

reference=$ref_dir/ref.fasta

gatk CreateSequenceDictionary -R $reference
samtools faidx $reference

array=($(find $ref_dir/Outputs/bams/ -name *.dedup.bam -printf "%f "))
id_arr=("${array[@]%.dedup.bam}")

#indv GVCFs
gatk HaplotypeCaller -R ${reference} --ploidy 1 -ERC GVCF \
    --I $ref_dir/Outputs/bams/${id_arr[$SLURM_ARRAY_TASK_ID]}.dedup.bam \
    -O $out_dir/indv_GVCFs/${id_arr[$SLURM_ARRAY_TASK_ID]}.raw.gvcf

#gatk HaplotypeCaller -R ${reference} --ploidy 1 -ERC GVCF \
#    --I $ref_dir/Outputs/bams/${id_arr[1]}.dedup.bam \
#    -O $out_dir/indv_GVCFs/${id_arr[1]}.raw.gvcf
