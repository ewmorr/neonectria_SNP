cd
mkdir core_fugr_invariant_sites_GVCF
mkdir core_fugr_invariant_sites_GVCF/indv_GVCFs

sbatch ~/repo/neonectria_SNP/pop_gen/phylogeny/indv_GVCFs_array.slurm

sbatch ~/repo/neonectria_SNP/pop_gen/phylogeny/genotype_gvcfs_invariant_sites.slurm

cd
mkdir core_neonectria_invariant_sites_GVCF
mkdir core_neonectria_invariant_sites_GVCF/indv_GVCFs

sbatch ~/repo/neonectria_SNP/pop_gen/phylogeny/indv_GVCFs_array.slurm

sbatch ~/repo/neonectria_SNP/pop_gen/phylogeny/genotype_gvcfs_invariant_sites.slurm

