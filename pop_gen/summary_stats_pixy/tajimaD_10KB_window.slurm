#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name="pixy"
#SBATCH --output=pixy.out
#SBATCH --partition=shared
#SBATCH --exclude=node105,node106,node107,node108,node117,node118,node119,node167

#mamba install danielecook::vcf-kit=0.2.6   "bwa>=0.7.17"   "samtools>=1.10"   "bcftools>=1.10"   "blast>=2.2.31"   "muscle>=3.8.31"   "primer3>=2.5.0" "gsl==2.5"
#need the last gsl call otherwise bcftools gets wrong library

module purge
module load anaconda/colsa
conda activate vcf-kit


cd ~/Nf_tajD/

arr=( $(cut -f 2 Nf.tajD_pops.tsv | sort | uniq) )

for i in ${arr[@]}
do 
    grep $i Nf.tajD_pops.tsv | cut -f 1 > tmp.sample_list
    bcftools view -S tmp.sample_list FINAL_snp.biallele.TajDdiploid.vcf.gz -Oz -o tmp.vcf.gz
    vk tajima 1E5 1E5 tmp.vcf.gz > $i.tajima.tsv
    rm tmp.vcf.gz
    rm tmp.sample_list
done

