conda activate bcftools

cd ~/repo/neonectria_SNP/data/Nf/final_tables/rm_dups/

#Convert GVCF to SNP matrix
gatk IndexFeatureFile --input FINAL_snp.IBD_analyses.vcf.gz
gatk VariantsToTable -V FINAL_snp.IBD_analyses.vcf.gz \
    -O FINAL_snp.IBD_analyses.table \
    -F CHROM -F POS -F REF -F ALT -F TYPE -GF GT 

#convert NA to N for fasta conversion
sed 's:\.:N:g' FINAL_snp.IBD_analyses.table > FINAL_snp.IBD_analyses.table.na2n

#there are INDELS and "MIXED" indel/snp sites
cut -f 5 FINAL_snp.IBD_analyses.table.na2n | sort | uniq

grep INDEL FINAL_snp.IBD_analyses.table.na2n | wc -l
# 107038
grep MIXED FINAL_snp.IBD_analyses.table.na2n | wc -l
# 12012
grep SNP FINAL_snp.IBD_analyses.table.na2n | wc -l
# 997091
#
#We will only deal with the SNPs for the phylogeny
grep -vE 'MIXED|INDEL' FINAL_snp.IBD_analyses.table.na2n > FINAL_snp.IBD_analyses.table.snps_only

# fasta conversion
perl ~/repo/neonectria_SNP/library/snp_table2fasta.pl FINAL_snp.IBD_analyses.table.snps_only 5 > FINAL_snp.snps_only.for_phylogeny.fasta
sed -i '' 's/NGT//' FINAL_snp.snps_only.for_phylogeny.fasta 

#also for invariants
#Convert GVCF to SNP matrix
gatk IndexFeatureFile --input FINAL_invariant.nuclear.vcf.gz
gatk VariantsToTable -V FINAL_invariant.nuclear.vcf.gz \
    -O FINAL_invariant.IBD_analyses.table \
    -F CHROM -F POS -F REF -F ALT -F TYPE -GF GT

#convert NA to N for fasta conversion
sed 's:\.:N:g' FINAL_invariant.IBD_analyses.table > FINAL_invariant.IBD_analyses.table.na2n

#there are INDELS and "MIXED" indel/snp sites
cut -f 5 FINAL_invariant.IBD_analyses.table.na2n | sort | uniq

grep INDEL FINAL_invariant.IBD_analyses.table.na2n | wc -l
# 106780
grep MIXED FINAL_invariant.IBD_analyses.table.na2n | wc -l
# 11811
grep SNP FINAL_invariant.IBD_analyses.table.na2n | wc -l
# 991996
#
#We will only deal with the SNPs for the phylogeny
grep -vE 'MIXED|INDEL' FINAL_invariant.IBD_analyses.table.na2n > FINAL_invariant.IBD_analyses.table.snps_only

# fasta conversion
perl ~/repo/neonectria_SNP/library/snp_table2fasta.pl FINAL_invariant.IBD_analyses.table.snps_only 5 > FINAL_invariant.snps_only.for_IBD.fasta
sed -i '' 's/NGT//' FINAL_invariant.snps_only.for_IBD.fasta

#Calculate ML tree in R
ml_tree.adegenet.r
ml_tree.adegenet.plot.r
#low div in VA and most related to WV then NC

cd ~/repo/neonectria_SNP/data/Nd/final_tables/rm_dups/

#Convert GVCF to SNP matrix
gatk IndexFeatureFile --input FINAL_snp.IBD_analyses.vcf.gz
gatk VariantsToTable -V FINAL_snp.IBD_analyses.vcf.gz \
    -O FINAL_snp.IBD_analyses.table \
    -F CHROM -F POS -F REF -F ALT -F TYPE -GF GT 

#convert NA to N for fasta conversion
sed 's:\.:N:g' FINAL_snp.IBD_analyses.table > FINAL_snp.IBD_analyses.table.na2n

#there are INDELS and "MIXED" indel/snp sites
cut -f 5 FINAL_snp.IBD_analyses.table.na2n | sort | uniq

grep INDEL FINAL_snp.IBD_analyses.table.na2n | wc -l
# 164793
grep MIXED FINAL_snp.IBD_analyses.table.na2n | wc -l
# 21564
grep SNP FINAL_snp.IBD_analyses.table.na2n | wc -l
# 1413299
#
#We will only deal with the SNPs for the phylogeny
grep -vE 'MIXED|INDEL' FINAL_snp.IBD_analyses.table.na2n > FINAL_snp.IBD_analyses.table.snps_only

# fasta conversion
perl ~/repo/neonectria_SNP/library/snp_table2fasta.pl FINAL_snp.IBD_analyses.table.snps_only 5 > FINAL_snp.snps_only.for_phylogeny.fasta
sed -i '' 's/NGT//' FINAL_snp.snps_only.for_phylogeny.fasta 

#also for invariants
#Convert GVCF to SNP matrix
gatk IndexFeatureFile --input FINAL_invariant.nuclear.vcf.gz
gatk VariantsToTable -V FINAL_invariant.nuclear.vcf.gz \
    -O FINAL_invariant.IBD_analyses.table \
    -F CHROM -F POS -F REF -F ALT -F TYPE -GF GT

#convert NA to N for fasta conversion
sed 's:\.:N:g' FINAL_invariant.IBD_analyses.table > FINAL_invariant.IBD_analyses.table.na2n

#there are INDELS and "MIXED" indel/snp sites
cut -f 5 FINAL_invariant.IBD_analyses.table.na2n | sort | uniq

grep INDEL FINAL_invariant.IBD_analyses.table.na2n | wc -l
# 164182
grep MIXED FINAL_invariant.IBD_analyses.table.na2n | wc -l
# 21263
grep SNP FINAL_invariant.IBD_analyses.table.na2n | wc -l
# 1407462
#
#We will only deal with the SNPs for the phylogeny
grep -vE 'MIXED|INDEL' FINAL_invariant.IBD_analyses.table.na2n > FINAL_invariant.IBD_analyses.table.snps_only

# fasta conversion
perl ~/repo/neonectria_SNP/library/snp_table2fasta.pl FINAL_invariant.IBD_analyses.table.snps_only 5 > FINAL_invariant.snps_only.for_IBD.fasta
sed -i '' 's/NGT//' FINAL_invariant.snps_only.for_IBD.fasta
