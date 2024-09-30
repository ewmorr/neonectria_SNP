# note, we originally ran this by performing nucmer alignment (fugr, nefa, nedi, neco) then parsing the core alignment fasta of fugr and mapping against that. The results were very poor (species intermixed on the tree). We deleted these results and are now performing the mapping against the whole fugr genome with the intention to then just split out the snps/invariant based on the core alignment
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

