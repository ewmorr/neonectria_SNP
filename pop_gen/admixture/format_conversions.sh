conda activate bcftools

#we want to run this on the filtered SNP file including all possible SNPs, i.e., the one we label for IBD analysis
cd ~/repo/neonectria_SNP/data/Nf/final_tables/rm_dups/

plink --vcf FINAL_snp.mac_ge2.biallele.LD.structure_analyses.vcf.gz --recode12 --out FINAL_snp.admixture --allow-extra-chr
sed 's/tig//g' FINAL_snp.admixture.map | sed 's/_pilon//g' > FINAL_snp.admixture.digit.map
mv FINAL_snp.admixture.digit.map FINAL_snp.admixture.map


cd ~/repo/neonectria_SNP/data/Nd/final_tables/rm_dups/
plink --vcf FINAL_snp.mac_ge2.biallele.LD.structure_analyses.vcf.gz --recode12 --out FINAL_snp.admixture --allow-extra-chr
sed 's/LDPL//g' FINAL_snp.admixture.map | sed 's/\.1//g' > FINAL_snp.admixture.digit.map
mv FINAL_snp.admixture.digit.map FINAL_snp.admixture.map


cd ~/repo/neonectria_SNP/data/Nc/final_tables/
plink --vcf FINAL_snp.mac_ge2.biallele.LD.structure_analyses.vcf.gz --recode12 --out FINAL_snp.admixture --allow-extra-chr
sed 's/WPDF//g' FINAL_snp.admixture.map | sed 's/\.1//g' > FINAL_snp.admixture.digit.map
mv FINAL_snp.admixture.digit.map FINAL_snp.admixture.map
