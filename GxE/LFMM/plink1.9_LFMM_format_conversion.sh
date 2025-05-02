conda activate bcftools

cd ~/repo/neonectria_SNP/data/Nf/final_tables/rm_dups/

# the full table
plink --vcf FINAL_snp.mac_ge2.biallele.gwas_analyses.vcf.gz --recode 01 --missing-genotype 9 --out FINAL_snp.mac_ge2.biallele.gwas_analyses.recode --allow-extra-chr

wc -l FINAL_snp.mac_ge2.biallele.gwas_analyses.recode.ped

cut -d " " -f 7- FINAL_snp.mac_ge2.biallele.gwas_analyses.recode.ped | awk  '{for (i=1;i<=NF;i+=2) printf "%s ", $i; printf "\n" }' > FINAL_snp.mac_ge2.biallele.gwas_analyses.lfmm

cut -d " " -f 1 FINAL_snp.mac_ge2.biallele.gwas_analyses.recode.ped > FINAL_snp.mac_ge2.biallele.gwas_analyses.sampleIDs

zgrep "##contig=<ID=" FINAL_snp.mac_ge2.biallele.gwas_analyses.vcf.gz > scaffold_lengths.txt

sed -e 's/##contig=<ID=//' -e 's/length=//' -e 's/>//' scaffold_lengths.txt > scaffold_lengths.csv

# the table including only cluster 2 (i.e., the homogeneous ancenstry cluster)

plink --vcf FINAL_snp.gwas_analyses.cluster2.vcf.gz --recode 01 --missing-genotype 9 --out FINAL_snp.gwas_analyses.cluster2.recode --allow-extra-chr

wc -l FINAL_snp.gwas_analyses.cluster2.recode.ped

cut -d " " -f 7- FINAL_snp.gwas_analyses.cluster2.recode.ped | awk  '{for (i=1;i<=NF;i+=2) printf "%s ", $i; printf "\n" }' > FINAL_snp.gwas_analyses.cluster2.lfmm

cut -d " " -f 1 FINAL_snp.gwas_analyses.cluster2.recode.ped > FINAL_snp.gwas_analyses.cluster2.sampleIDs


# the table including only cluster 2 (i.e., the homogeneous ancenstry cluster)
# with no NC samples

plink --vcf FINAL_snp.gwas_analyses.cluster2_no_NC.vcf.gz --recode 01 --missing-genotype 9 --out FINAL_snp.gwas_analyses.cluster2_no_NC.recode --allow-extra-chr

wc -l FINAL_snp.gwas_analyses.cluster2_no_NC.recode.ped

cut -d " " -f 7- FINAL_snp.gwas_analyses.cluster2_no_NC.recode.ped | awk  '{for (i=1;i<=NF;i+=2) printf "%s ", $i; printf "\n" }' > FINAL_snp.gwas_analyses.cluster2_no_NC.lfmm

cut -d " " -f 1 FINAL_snp.gwas_analyses.cluster2_no_NC.recode.ped > FINAL_snp.gwas_analyses.cluster2_no_NC.sampleIDs

