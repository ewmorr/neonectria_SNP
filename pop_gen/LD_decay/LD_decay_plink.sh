#trying this locally
# run on four cores (we have 5 performance and 6 efficiency cores)

conda activate bcftools

#we want to run this on the filtered SNP file including all possible SNPs, i.e., the one we label for IBD analysis
cd ~/repo/neonectria_SNP/data/Nf/final_tables/rm_dups/

#bfile takes bed, bim, fam files
#we already converted this for calculation of hamming dist
# --ld-window 10-1 --ld-window-kb 1000 #these are the defaults
# --ld-window seems to control the number of pairwise comparisons examined/output NOT the minimum gen distance
# --threads 8
plink --bfile FINAL_snp.IBD_analyses --r2 gz --ld-window-r2 0 --threads 8 --allow-extra-chr --out FINAL_snp.LD_decay.retry
# plink --bfile FINAL_snp.IBD_analyses --r2 gz --ld-window-r2 0 --ld-window 999999 --ld-window-kb 6000 --threads 8 --allow-extra-chr --out FINAL_snp.LD_decay #these files are HUGE! like 100-300 Gb. The 10K comparison version is smaller, 20-30 Gb. We can keep that v for now but let's try with 1K comparisons first
plink --bfile FINAL_snp.IBD_analyses --r2 gz --ld-window-r2 0 --ld-window 9999 --ld-window-kb 6000 --threads 8 --allow-extra-chr --out FINAL_snp.LD_decay_9999
plink --bfile FINAL_snp.IBD_analyses --r2 gz --ld-window-r2 0 --ld-window 999 --ld-window-kb 6000 --threads 8 --allow-extra-chr --out FINAL_snp.LD_decay_999 # this is about 4 Gb. Much more manageable
#1116141 variants and 115 people pass filters and QC.

cd ~/repo/neonectria_SNP/data/Nd/final_tables/rm_dups/
plink --bfile FINAL_snp.IBD_analyses --r2 gz --ld-window-r2 0 --threads 8 --allow-extra-chr --out FINAL_snp.LD_decay.retry
#plink --bfile FINAL_snp.IBD_analyses --r2 gz --ld-window-r2 0 --ld-window 999999 --ld-window-kb 6000 --threads 8 --allow-extra-chr --out FINAL_snp.LD_decay
plink --bfile FINAL_snp.IBD_analyses --r2 gz --ld-window-r2 0 --ld-window 9999 --ld-window-kb 6000 --threads 8 --allow-extra-chr --out FINAL_snp.LD_decay_9999
plink --bfile FINAL_snp.IBD_analyses --r2 gz --ld-window-r2 0 --ld-window 999 --ld-window-kb 6000 --threads 8 --allow-extra-chr --out FINAL_snp.LD_decay_999
#1599656 variants and 30 people pass filters and QC.

cd ~/repo/neonectria_SNP/data/Nc/final_tables/
plink --bfile FINAL_snp.IBD_analyses --r2 gz --ld-window-r2 0 --threads 8 --allow-extra-chr --out FINAL_snp.LD_decay.retry
#plink --bfile FINAL_snp.IBD_analyses --r2 gz --ld-window-r2 0 --ld-window 999999 --ld-window-kb 6000 --threads 8 --allow-extra-chr --out FINAL_snp.LD_decay
plink --bfile FINAL_snp.IBD_analyses --r2 gz --ld-window-r2 0 --ld-window 9999 --ld-window-kb 6000 --threads 8 --allow-extra-chr --out FINAL_snp.LD_decay_9999
plink --bfile FINAL_snp.IBD_analyses --r2 gz --ld-window-r2 0 --ld-window 999 --ld-window-kb 6000 --threads 8 --allow-extra-chr --out FINAL_snp.LD_decay_999
#416775 variants and 5 people pass filters and QC.
