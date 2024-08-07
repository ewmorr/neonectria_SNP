conda activate bcftools

cd ~/repo/neonectria_SNP/data/Nf/final_tables/rm_dups/

#Convert GVCF to SNP matrix
gatk IndexFeatureFile --input FINAL_snp.IBD_analyses.vcf.gz
gatk VariantsToTable -V FINAL_snp.IBD_analyses.vcf.gz \
    -O FINAL_snp.IBD_analyses.table \
    -F CHROM -F POS -F REF -F ALT -F TYPE -GF GT 

head -n 1 FINAL_snp.IBD_analyses.table | grep -o "NG" | wc -l
#convert NA to N for fasta conversion
sed 's:\.:N:g' FINAL_snp.IBD_analyses.table > FINAL_snp.IBD_analyses.table.na2n
head -n 1 FINAL_snp.IBD_analyses.table.na2n | grep -o "NGT" | wc -l

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
head -n 1 FINAL_snp.IBD_analyses.table.snps_only | grep -o "NGT" | wc -l

# fasta conversion
perl ~/repo/neonectria_SNP/library/snp_table2fasta.pl FINAL_snp.IBD_analyses.table.snps_only 5 > FINAL_snp.snps_only.for_phylogeny.fasta
sed -i '' 's/NGT//' FINAL_snp.snps_only.for_phylogeny.fasta 

#also for invariants
#Convert GVCF to SNP matrix
gatk IndexFeatureFile --input FINAL_invariant.nuclear.vcf.gz
gatk VariantsToTable -V FINAL_invariant.nuclear.vcf.gz \
    -O FINAL_invariant.IBD_analyses.table \
    -F CHROM -F POS -F REF -F ALT -F TYPE -GF GT
head -n 1 FINAL_invariant.IBD_analyses.table | grep -o "NG" | wc -l
#115 

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

#We need to split the table because ape::dist.dna maxes at about 2.3*10^9 bases
# split by contig
#note we tested and either gaps or Ns are ignored in a pairwise fashion by dist.dna

cut -f 1 FINAL_invariant.IBD_analyses.table.snps_only | sort | uniq > tigs_names.txt
mkdir invariant_table_tigs

while IFS= read -r line
do
    grep $line FINAL_invariant.IBD_analyses.table.snps_only > invariant_table_tigs/$line.table.snps_only &
done < tigs_names.txt

cd invariant_table_tigs
for i in tig*
do(
    #echo $i
    #cat CHROM.table.snps_only $i > $i.tmp & mv $i.tmp $i
    wc -l $i
)
done

#the lens seem odd so retrying
mkdir invariant_table_tigs_retry

while IFS= read -r line
do
    grep $line FINAL_invariant.IBD_analyses.table.snps_only > invariant_table_tigs_retry/$line.table.snps_only &
done < tig_names_1.txt

while IFS= read -r line
do
    grep $line FINAL_invariant.IBD_analyses.table.snps_only > invariant_table_tigs_retry/$line.table.snps_only &
done < tig_names_2.txt

 grep tig00000320_pilon FINAL_invariant.IBD_analyses.table.snps_only > invariant_table_tigs_retry/tig00000320_pilon.table.snps_only 

cd invariant_table_tigs_retry
rm lens.txt
for i in tig*only
do(
    #echo $i
    #cat CHROM.table.snps_only $i > $i.tmp 
    mv $i.tmp $i
    wc -l $i >> lens_cat.txt
)
done
# seems like the first round got messed up, maybe from overwriting or the processes crashing, it is about 6M bp too big
# This round is the correct len

# cd ..
# rm -r invariant_table_tigs
# mv invariant_table_tigs_retry invariant_table_tigs
# 
# 
# run each tig on server `~/repo/neonectria_SNP/premise/invariant_hamming.Nf.slurm`
# doesn't need extra high mem when run tig wise but should grab 100G to be safe (pulled in 40-50 when wewathced the first tig)

# fasta conversion
# the following needs to be run on a high mem. Took 350G of memory on initial run. Should fix the script
#perl ~/repo/neonectria_SNP/library/snp_table2fasta.pl FINAL_invariant.IBD_analyses.table.snps_only 5 > FINAL_invariant.snps_only.for_IBD.fasta
#sed -i '' 's/NGT//' FINAL_invariant.snps_only.for_IBD.fasta

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
