conda activate plink1.9

#Nf
cd ~/repo/neonectria_SNP/data/Nf/final_tables/rm_dups/
gunzip FINAL_snp.IBD_analyses.vcf.gz 

plink --vcf FINAL_snp.IBD_analyses.vcf --out FINAL_snp.IBD_analyses --allow-extra-chr --make-bed
bgzip FINAL_snp.IBD_analyses.vcf
plink --bfile FINAL_snp.IBD_analyses --allow-extra-chr --distance 
mv plink.dist FINAL_snp.IBD_analyses.dist
mv plink.dist.id FINAL_snp.IBD_analyses.dist.id

bcftools view FINAL_invariant.nuclear.vcf.gz -Ob -o FINAL_invariant.nuclear.bcf
plink --bcf FINAL_invariant.nuclear.bcf --out FINAL_invariant.IBD_analyses --allow-extra-chr --make-bed --threads 8
plink --bfile FINAL_invariant.IBD_analyses --allow-extra-chr --distance square --threads 8
mv plink.dist FINAL_invariant.IBD_analyses.dist
mv plink.dist.id FINAL_invariant.IBD_analyses.dist.id

#Nd
cd ~/repo/neonectria_SNP/data/Nd/final_tables/rm_dups/

bcftools view FINAL_snp.IBD_analyses.vcf.gz -Ob -o FINAL_snp.IBD_analyses.bcf
plink --bcf FINAL_snp.IBD_analyses.bcf --out FINAL_snp.IBD_analyses --allow-extra-chr --make-bed
plink --bfile FINAL_snp.IBD_analyses --allow-extra-chr --distance 
mv plink.dist FINAL_snp.IBD_analyses.dist
mv plink.dist.id FINAL_snp.IBD_analyses.dist.id

bcftools view FINAL_invariant.nuclear.vcf.gz -Ob -o FINAL_invariant.nuclear.bcf
plink --bcf FINAL_invariant.nuclear.bcf --out FINAL_invariant.IBD_analyses --allow-extra-chr --make-bed --threads 8
plink --bfile FINAL_invariant.IBD_analyses --allow-extra-chr --distance square --threads 8
mv plink.dist FINAL_invariant.IBD_analyses.dist
mv plink.dist.id FINAL_invariant.IBD_analyses.dist.id

#Nc
cd ~/repo/neonectria_SNP/data/Nc/final_tables/

bcftools view FINAL_snp.IBD_analyses.vcf.gz -Ob -o FINAL_snp.IBD_analyses.bcf
plink --bcf FINAL_snp.IBD_analyses.bcf --out FINAL_snp.IBD_analyses --allow-extra-chr --make-bed
plink --bfile FINAL_snp.IBD_analyses --allow-extra-chr --distance 
mv plink.dist FINAL_snp.IBD_analyses.dist
mv plink.dist.id FINAL_snp.IBD_analyses.dist.id

bcftools view FINAL_invariant.nuclear.vcf.gz -Ob -o FINAL_invariant.nuclear.bcf
plink --bcf FINAL_invariant.nuclear.bcf --out FINAL_invariant.IBD_analyses --allow-extra-chr --make-bed --threads 8
plink --bfile FINAL_invariant.IBD_analyses --allow-extra-chr --distance square --threads 8
mv plink.dist FINAL_invariant.IBD_analyses.dist
mv plink.dist.id FINAL_invariant.IBD_analyses.dist.id
