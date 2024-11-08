#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J bcftools
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output bcftools_core_filter.out
#SBATCH --exclude=node117,node118

module purge
module load linuxbrew/colsa

cd $HOME/core_genome_fugr_SPANDx/Outputs/Master_vcf/

bgzip -c out.vcf > out.vcf.gz
tabix out.vcf.gz

bcftools view -R $HOME/core_fugr_invariant_sites_GVCF/vcf.collated.coords -Oz -o core.out.vcf.gz out.vcf.gz
