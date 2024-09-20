#!/usr/bin/bash
#SBATCH --job-name="gatk"
#SBATCH --output=gatk.out
#SBATCH --partition=shared
#SBATCH --ntasks=1 --cpus-per-task=1

module purge
module load anaconda/colsa

conda activate spandx

#First call indv_GVCFs_array.slurm for individual vcf files

out_dir=/mnt/home/garnas/ewj4/Nc_invariant_sites_GVCF
in_dir=/mnt/home/garnas/ewj4/Nc_invariant_sites_GVCF/indv_GVCFs

reference=/mnt/home/garnas/ewj4/Nc_SPANDx_all_seqs/ref.fasta

#SPANDx gatk configs
    CLUSTER_SNP=3
    CLUSTER_WINDOW_SNP=10
    MLEAF_SNP=0.95
    QD_SNP=10.0
    MQ_SNP=30.0
    FS_SNP=60.0
    QUAL_SNP=30.0
    LOW_DEPTH=2
    HIGH_DEPTH=3
#end

array=($(find $in_dir/ -name *.gvcf -printf "%f "))
array2=("${array[@]/#/-V $in_dir/}")

gatk CombineGVCFs -R $reference ${array2[*]} \
    -O $out_dir/master.invariant_sites.vcf
gatk GenotypeGVCFs -R $reference \
    -V $out_dir/master.invariant_sites.vcf \
    -all-sites \
    -O $out_dir/out.invariant_sites.vcf

gatk VariantFiltration -R ${reference} \
    -O $out_dir/out.invariant_sites.filtered.vcf \
    -V $out_dir/out.invariant_sites.vcf \
    --cluster-size $CLUSTER_SNP -window $CLUSTER_WINDOW_SNP \
    -filter "QD < $QD_SNP" --filter-name "QDFilter" \
    -filter "MQ < $MQ_SNP" --filter-name "MQFilter" \
    -filter "FS > $FS_SNP" --filter-name "HaplotypeScoreFilter"
