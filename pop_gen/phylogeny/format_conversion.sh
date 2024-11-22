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


##############################
#Nc

cd ~/repo/neonectria_SNP/data/Nc/final_tables/

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
# 36883
grep MIXED FINAL_snp.IBD_analyses.table.na2n | wc -l
# 1109
grep SNP FINAL_snp.IBD_analyses.table.na2n | wc -l
# 378783
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
# 36883
grep MIXED FINAL_invariant.IBD_analyses.table.na2n | wc -l
# 1109
grep SNP FINAL_invariant.IBD_analyses.table.na2n | wc -l
# 378785
#
#We will only deal with the SNPs for the phylogeny
grep -vE 'MIXED|INDEL' FINAL_invariant.IBD_analyses.table.na2n > FINAL_invariant.IBD_analyses.table.snps_only

# fasta conversion
#perl ~/repo/neonectria_SNP/library/snp_table2fasta.pl FINAL_invariant.IBD_analyses.table.snps_only 5 > FINAL_invariant.snps_only.for_IBD.fasta
#sed -i '' 's/NGT//' FINAL_invariant.snps_only.for_IBD.fasta


##############################
#shared_buscos

cd ~/repo/neonectria_SNP/data/shared_buscos/final_tables/rm_dups

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
# 42022
grep MIXED FINAL_snp.IBD_analyses.table.na2n | wc -l
# 4722
grep SNP FINAL_snp.IBD_analyses.table.na2n | wc -l
# 695814
#
#We will only deal with the SNPs for the phylogeny
grep -vE 'MIXED|INDEL' FINAL_snp.IBD_analyses.table.na2n > FINAL_snp.IBD_analyses.table.snps_only

# fasta conversion
perl ~/repo/neonectria_SNP/library/snp_table2fasta.pl FINAL_snp.IBD_analyses.table.snps_only 5 > FINAL_snp.snps_only.for_phylogeny.fasta
sed -i '' 's/NGT//' FINAL_snp.snps_only.for_phylogeny.fasta # would need to fix this to make sure we are not removing NGT sequences but we are now using the core genome method

##############################
#core_fugr / Fugr1_ref

cd ~/repo/neonectria_SNP/data/Fugr1_ref/final_tables

#Convert GVCF to SNP matrix
bcftools view -Oz -o FINAL_invariant.vcf.gz FINAL_invariant.bcf
gatk IndexFeatureFile --input FINAL_invariant.vcf.gz
gatk VariantsToTable -V FINAL_invariant.vcf.gz \
    -O FINAL_invariant.IBD_analyses.table \
    -F CHROM -F POS -F REF -F ALT -F TYPE -GF GT

#convert NA to N for fasta conversion
sed 's:\t\.:\tN:g' FINAL_invariant.IBD_analyses.table > FINAL_invariant.IBD_analyses.table.na2n

#there are INDELS and "MIXED" indel/snp sites
cut -f 5 FINAL_invariant.IBD_analyses.table.na2n | sort | uniq

grep INDEL FINAL_invariant.IBD_analyses.table.na2n | wc -l
# 5034
grep MIXED FINAL_invariant.IBD_analyses.table.na2n | wc -l
# 1276
grep SNP FINAL_invariant.IBD_analyses.table.na2n | wc -l
# 35792
#
#We will only deal with the SNPs and INVARIANT for the phylogeny
grep -vE 'MIXED|INDEL' FINAL_invariant.IBD_analyses.table.na2n > FINAL_invariant.IBD_analyses.table.snps_only

# fasta conversion
perl ~/repo/neonectria_SNP/library/snp_table2fasta.pl FINAL_invariant.IBD_analyses.table.snps_only 5 > FINAL_invariant.snps_only.for_phylogeny.fasta
#sed -i '' 's/NGT//' FINAL_invariant.snps_only.for_phylogeny.fasta # note if this is needed need to correct the regex to only look in the header

#############################
#############################
### Need to add FUGR ref seq to the fasta manually
### Pull the list of positions from the final table (that was convertedc to fasta)
### 
### 
grep -v "^#" FINAL_invariant.IBD_analyses.table.na2n | cut -f 1,2 | less
grep -v "^#" FINAL_invariant.IBD_analyses.table.snps_only | cut -f 1,2 | less
grep -v "^#" FINAL_invariant.IBD_analyses.table.snps_only | cut -f 1,2 > FINAL_invariant.IBD_analyses.table.snps_only.positions

cd ~/repo/neonectria_SNP/data/Fugr1_ref

perl ../../library/get_fasta_from_pos.pl core.Fusgr1_AssemblyScaffolds_Repeatmasked-N.fasta final_tables/FINAL_invariant.IBD_analyses.table.snps_only.positions > core.Fusgr1.snp_pos.fasta

#check to make sure the tig orders are correct before conactenate
cut -f 1 final_tables/FINAL_invariant.IBD_analyses.table.snps_only.positions | uniq | sed '1d' > pos_scf_order.txt
grep ">" core.Fusgr1.snp_pos.fasta | sed 's/^>//' > fas_scf_order.txt

wc -l pos_scf_order.txt
wc -l fas_scf_order.txt
#both 430

diff pos_scf_order.txt fas_scf_order.txt
# no diffs

perl ../../library/cat_fasta.pl core.Fusgr1.snp_pos.fasta Fusgr1_core_snp_pos > core.Fusgr1.snp_pos.cat.fasta
perl ../../library/get_seq_lens.pl core.Fusgr1.snp_pos.cat.fasta
#290590
perl ../../library/get_seq_lens.pl final_tables/FINAL_invariant.snps_only.for_phylogeny.fasta | less
#290590

sed '1d' final_tables/FINAL_invariant.IBD_analyses.table.snps_only.positions | wc -l
# 290590

cat core.Fusgr1.snp_pos.cat.fasta final_tables/FINAL_invariant.snps_only.for_phylogeny.fasta > core.Fusgr1-neonectria.snps_plus_invariant_aln.fasta
perl ../../library/get_seq_lens.pl core.Fusgr1-neonectria.snps_plus_invariant_aln.fasta | less
# looks good



##############################
# no core extract/ Fugr1_ref

cd ~/repo/neonectria_SNP/data/Fugr1_ref/map_against_Fugr_no_core_extract/final_tables

#Convert GVCF to SNP matrix
bcftools view -Oz -o FINAL_snp.vcf.gz FINAL_snp.bcf
gatk IndexFeatureFile --input FINAL_snp.vcf.gz
gatk VariantsToTable -V FINAL_snp.vcf.gz \
    -O FINAL_snp.IBD_analyses.table \
    -F CHROM -F POS -F REF -F ALT -F TYPE -GF GT

#convert NA to N for fasta conversion
sed 's:\t\.:\tN:g' FINAL_snp.IBD_analyses.table > FINAL_snp.IBD_analyses.table.na2n

#there are INDELS and "MIXED" indel/snp sites
cut -f 5 FINAL_snp.IBD_analyses.table.na2n | sort | uniq

grep INDEL FINAL_snp.IBD_analyses.table.na2n | wc -l
# 80841
grep MIXED FINAL_snp.IBD_analyses.table.na2n | wc -l
# 19042
grep SNP FINAL_snp.IBD_analyses.table.na2n | wc -l
# 290104 (compared to 34K from the core aln)
#
#We will only deal with the SNPs for the phylogeny
grep -vE 'MIXED|INDEL' FINAL_snp.IBD_analyses.table.na2n > FINAL_snp.IBD_analyses.table.snps_only

# fasta conversion
perl ~/repo/neonectria_SNP/library/snp_table2fasta.pl FINAL_snp.IBD_analyses.table.snps_only 5 > FINAL_snp.snps_only.for_phylogeny.fasta
less -S FINAL_snp.snps_only.for_phylogeny.fasta
#sed -i '' 's/NGT//' FINAL_snp.snps_only.for_phylogeny.fasta # note if this is needed need to correct the regex to only look in the header

#############################
#############################
### Need to add FUGR ref seq to the fasta manually
### Pull the list of positions from the final table (that was converted to fasta)
###
###
grep -v "^#" FINAL_snp.IBD_analyses.table.na2n | cut -f 1,2 | less
grep -v "^#" FINAL_snp.IBD_analyses.table.snps_only | cut -f 1,2 | less
grep -v "^#" FINAL_snp.IBD_analyses.table.snps_only | cut -f 1,2 > FINAL_snp.IBD_analyses.table.snps_only.positions

cd ~/repo/neonectria_SNP/data/Fugr1_ref

#note this MUST be against the full ref, not the nucmer core ref (otherwise the positions will not be calculated correctly in the VCF)
#i.e., the ref must be the same as what was mapped against to create the VCF
perl ../../library/get_fasta_from_pos.pl Fusgr1_AssemblyScaffolds_Repeatmasked-N.fasta map_against_Fugr_no_core_extract/final_tables/FINAL_snp.IBD_analyses.table.snps_only.positions > map_against_Fugr_no_core_extract/core.Fusgr1.snp_pos.fasta

cd map_against_Fugr_no_core_extract
#check to make sure the tig orders are correct before conactenate
cut -f 1 final_tables/FINAL_snp.IBD_analyses.table.snps_only.positions | uniq | sed '1d' > pos_scf_order.txt
grep ">" core.Fusgr1.snp_pos.fasta | sed 's/^>//' > fas_scf_order.txt

wc -l pos_scf_order.txt
wc -l fas_scf_order.txt
#both 15

diff pos_scf_order.txt fas_scf_order.txt
# no diffs

perl ../../../library/cat_fasta.pl core.Fusgr1.snp_pos.fasta Fusgr1_core_snp_pos > core.Fusgr1.snp_pos.cat.fasta
perl ../../../library/get_seq_lens.pl core.Fusgr1.snp_pos.cat.fasta
#290104 (funny enough just a bit shorter than the core SNPs plus invariants)
perl ../../../library/get_seq_lens.pl final_tables/FINAL_snp.snps_only.for_phylogeny.fasta | less
#290104

sed '1d' final_tables/FINAL_snp.IBD_analyses.table.snps_only.positions | wc -l
# 290104

cat core.Fusgr1.snp_pos.cat.fasta final_tables/FINAL_snp.snps_only.for_phylogeny.fasta > core.Fusgr1-neonectria.snps_aln.fasta
perl ../../../library/get_seq_lens.pl core.Fusgr1-neonectria.snps_aln.fasta | less
# looks good
