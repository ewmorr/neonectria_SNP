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

# let's try with thinning right off the bat
plink --bfile FINAL_snp.IBD_analyses --r2 gz --ld-window-r2 0 --thin 0.01 --ld-window 9999 --ld-window-kb 6000 --threads 8 --allow-extra-chr --out FINAL_snp.LD_decay_9999
#11225 variants (1%)
gunzip -c FINAL_snp.LD_decay_9999.ld.gz | tail -n +2 | tr -s ' ' | cut -d ' ' -f 2,5,7 > FINAL_snp.LD_decay_9999.ld.apos-bpos-r
wc -l FINAL_snp.LD_decay_9999.ld.apos-bpos-r
#3330462 very manageable

#after plotting the above the loess fit minimizes quickly way below 1Mb. We will restrict window size to 1000
plink --bfile FINAL_snp.IBD_analyses --r2 gz --ld-window-r2 0 --thin 0.01 --ld-window 9999 --ld-window-kb 1000 --threads 8 --allow-extra-chr --out FINAL_snp.1Mb_max.LD_decay_9999
gunzip -c FINAL_snp.1Mb_max.LD_decay_9999.ld.gz | tail -n +2 | tr -s ' ' | cut -d ' ' -f 2,5,7 > FINAL_snp.1Mb_max.LD_decay_9999.ld.apos-bpos-r
wc -l FINAL_snp.1Mb_max.LD_decay_9999.ld.apos-bpos-r
#11078 variants
#1915391

#plink --bfile FINAL_snp.IBD_analyses --r2 gz --ld-window-r2 0 --ld-window 999 --ld-window-kb 6000 --threads 8 --allow-extra-chr --out FINAL_snp.LD_decay_999 # this is about 4 Gb. Much more manageable
#1116141 variants and 115 people pass filters and QC.

#gunzip -c FINAL_snp.LD_decay_999.ld.gz | tail -n +2 | tr -s ' ' | cut -d ' ' -f 2,5,7 > FINAL_snp.LD_decay_999.ld.apos-bpos-r
# the full file is 14.5 G. We will gzip this but also random sample down to maybe 1M lines. Presumably the full file is somewhere around 1.1M*999 lines
# the 1M random sample is 26 Mb. We can increase this to maybe 2.5 M

#shuf -n 2500000 FINAL_snp.LD_decay_999.ld.apos-bpos-r > FINAL_snp.LD_decay_999.ld.apos-bpos-r.2.5M &
# the 2.5 M is 65 Mb. Let's go with 10 M

#shuf -n 10000000 FINAL_snp.LD_decay_999.ld.apos-bpos-r > FINAL_snp.LD_decay_999.ld.apos-bpos-r.10M &
#gzip FINAL_snp.LD_decay_999.ld.apos-bpos-r && rm FINAL_snp.LD_decay_999.ld.apos-bpos-r


cd ~/repo/neonectria_SNP/data/Nd/final_tables/rm_dups/
plink --bfile FINAL_snp.IBD_analyses --r2 gz --ld-window-r2 0 --thin 0.01 --ld-window 9999 --ld-window-kb 6000 --threads 8 --allow-extra-chr --out FINAL_snp.LD_decay_9999
#15778 variants (1%)
gunzip -c FINAL_snp.LD_decay_9999.ld.gz | tail -n +2 | tr -s ' ' | cut -d ' ' -f 2,5,7 > FINAL_snp.LD_decay_9999.ld.apos-bpos-r
wc -l FINAL_snp.LD_decay_9999.ld.apos-bpos-r
#2810354

plink --bfile FINAL_snp.IBD_analyses --r2 gz --ld-window-r2 0 --thin 0.01 --ld-window 9999 --ld-window-kb 1000 --threads 8 --allow-extra-chr --out FINAL_snp.1Mb_max.LD_decay_9999
gunzip -c FINAL_snp.1Mb_max.LD_decay_9999.ld.gz | tail -n +2 | tr -s ' ' | cut -d ' ' -f 2,5,7 > FINAL_snp.1Mb_max.LD_decay_9999.ld.apos-bpos-r
wc -l FINAL_snp.1Mb_max.LD_decay_9999.ld.apos-bpos-r
#16016 variants
#1893111

cd ~/repo/neonectria_SNP/data/Nc/final_tables/
plink --bfile FINAL_snp.IBD_analyses --r2 gz --ld-window-r2 0 --thin 0.1 --ld-window 9999 --ld-window-kb 6000 --threads 8 --allow-extra-chr --out FINAL_snp.LD_decay_9999
#41559 variants (10%)
gunzip -c FINAL_snp.LD_decay_9999.ld.gz | tail -n +2 | tr -s ' ' | cut -d ' ' -f 2,5,7 > FINAL_snp.LD_decay_9999.ld.apos-bpos-r
wc -l FINAL_snp.LD_decay_9999.ld.apos-bpos-r
#3221473

plink --bfile FINAL_snp.IBD_analyses --r2 gz --ld-window-r2 0 --thin 0.05 --ld-window 9999 --ld-window-kb 1000 --threads 8 --allow-extra-chr --out FINAL_snp.1Mb_max.LD_decay_9999
gunzip -c FINAL_snp.1Mb_max.LD_decay_9999.ld.gz | tail -n +2 | tr -s ' ' | cut -d ' ' -f 2,5,7 > FINAL_snp.1Mb_max.LD_decay_9999.ld.apos-bpos-r
wc -l FINAL_snp.1Mb_max.LD_decay_9999.ld.apos-bpos-r
#20909 variants
#817469
