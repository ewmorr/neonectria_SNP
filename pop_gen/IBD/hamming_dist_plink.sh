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

bcftools view FINAL_snp.mac_ge2.LD.pca_analyses.vcf.gz -Ob -o FINAL_snp.mac_ge2.LD.pca_analyses.bcf
plink --bcf FINAL_snp.mac_ge2.LD.pca_analyses.bcf --out FINAL_snp.mac_ge2.LD.pca_analyses --allow-extra-chr --make-bed --threads 8
plink --bfile FINAL_snp.mac_ge2.LD.pca_analyses --allow-extra-chr --distance square --threads 8
mv plink.dist FINAL_snp.mac_ge2.LD.pca_analyses.dist
mv plink.dist.id FINAL_snp.mac_ge2.LD.pca_analyses.dist.id

gunzip -c FINAL_snp.mac_ge2.LD.pca_analyses.vcf.gz | grep -v ^# | wc -l
# 196928
gunzip -c FINAL_snp.mac_ge2.biallele.LD.structure_analyses.vcf.gz | grep -v ^# | wc -l
# 188237


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

bcftools view FINAL_snp.mac_ge2.LD.pca_analyses.vcf.gz -Ob -o FINAL_snp.mac_ge2.LD.pca_analyses.bcf
plink --bcf FINAL_snp.mac_ge2.LD.pca_analyses.bcf --out FINAL_snp.mac_ge2.LD.pca_analyses --allow-extra-chr --make-bed --threads 8
plink --bfile FINAL_snp.mac_ge2.LD.pca_analyses --allow-extra-chr --distance square --threads 8
mv plink.dist FINAL_snp.mac_ge2.LD.pca_analyses.dist
mv plink.dist.id FINAL_snp.mac_ge2.LD.pca_analyses.dist.id

gunzip -c FINAL_snp.mac_ge2.LD.pca_analyses.vcf.gz | grep -v ^# | wc -l
# 87437
gunzip -c FINAL_snp.mac_ge2.biallele.LD.structure_analyses.vcf.gz | grep -v ^# | wc -l
# 81334

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

bcftools view FINAL_snp.mac_ge2.LD.pca_analyses.vcf.gz -Ob -o FINAL_snp.mac_ge2.LD.pca_analyses.bcf
plink --bcf FINAL_snp.mac_ge2.LD.pca_analyses.bcf --out FINAL_snp.mac_ge2.LD.pca_analyses --allow-extra-chr --make-bed --threads 8
plink --bfile FINAL_snp.mac_ge2.LD.pca_analyses --allow-extra-chr --distance square --threads 8
mv plink.dist FINAL_snp.mac_ge2.LD.pca_analyses.dist
mv plink.dist.id FINAL_snp.mac_ge2.LD.pca_analyses.dist.id

