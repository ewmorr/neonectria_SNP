#!/usr/bin/bash

#mamba install danielecook::vcf-kit=0.2.6   "bwa>=0.7.17"   "samtools>=1.10"   "bcftools>=1.10"   "blast>=2.2.31"   "muscle>=3.8.31"   "primer3>=2.5.0" "gsl==2.5"
#need the last gsl call otherwise bcftools gets wrong library

# the below works on the head node on Premise but throws an error about a missing python module when run via slurm... weird

conda activate vcf-kit

cd ~/repo/neonectria_SNP/data/Nf/pixy

pops_file=~/repo/neonectria_SNP/data/sample_metadata/Nf.tajD_pops.tsv
vcf_file=~/repo/neonectria_SNP/data/Nf/final_tables/rm_dups/FINAL_snp.biallele.TajDdiploid.vcf.gz

arr=( $(cut -f 2 $pops_file | sort | uniq) )

for i in ${arr[@]}
do 
    grep $i $pops_file | cut -f 1 > tmp.sample_list
    bcftools view -S tmp.sample_list $vcf_file -Oz -o tmp.vcf.gz
    vk tajima 1E5 1E5 tmp.vcf.gz > $i.tajima.tsv
    rm tmp.vcf.gz
    rm tmp.sample_list
done

