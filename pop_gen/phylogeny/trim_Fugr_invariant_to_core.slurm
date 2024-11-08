#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J bcftools
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output bcftools_core_filter.out
#SBATCH --exclude=node117,node118

module purge
module load linuxbrew/colsa

cd $HOME/core_fugr_invariant_sites_GVCF

#bgzip -c out.invariant_sites.vcf > out.invariant_sites.vcf.gz
#tabix out.invariant_sites.vcf.gz

bcftools view -R vcf.collated.coords -Oz -o core.out.invariant_sites.vcf.gz out.invariant_sites.vcf.gz
