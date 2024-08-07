
#run locally
conda activate bcftools

cd ~/repo/neonectria_SNP/data/Nc/

QDparam=2
MQparam=40
FSparam=60

gatk CreateSequenceDictionary -R ref.fasta
samtools faidx ref.fasta

gatk VariantFiltration -R ref.fasta -V out.vcf -O INFOfilters.vcf \
    --filter-name "MQ40" -filter "MQ < 40.0" \
    --filter-name "QD2" -filter "QD < 2.0" \
    --filter-name "FS60" -filter "FS > 60.0"
    

#remove commas from vcf INFO field descriptions before proceeding with vcftools
# https://unix.stackexchange.com/questions/48672/only-remove-commas-embedded-within-quotes-in-a-comma-delimited-file

perl -pe 's/(".+?[^\\]")/($ret = $1) =~ (s#,##g); $ret/ge' INFOfilters.vcf > INFOfilters.noComma.vcf


# For the diversity analyses we will keep all mac and all allele counts. We have other qual fiulters to depend on to filter out bad calls, and these types of calls are obviously important to calculation of pi and Tajima's D, for example. Also see Lou et al 2021 which points out these should not be removed for div analyses

vcftools --vcf INFOfilters.noComma.vcf --out INFOfilters.removed --recode --remove-filtered-all
#kept 5 of 5 individuals
# After filtering, kept 488725 out of a possible 575775 Sites
# 575775 - 488725 = 87050

#There are no extreme %NA outliers so we just proceed with filtering as normal (as opposed to Nf)

####################
#DP filters
# we are prforming filtering on the INFO filed filtered dataset, but NOT having removed biallelic and MAC filters.
# This dataset will be used for diversity metrics. These are presumably good quality calls based on the other filters
#For any GWAS or GxE and STRUCTURE we will want to filter for biallele and MAC as the methods either ignore or
# don't handle unfiltered
# can perform biallelic and MAC  filters on the table after performing the DP and NA filters

# we ran DP filters on each library separately for Nf and Nd, but for Nc we only have the one, so we can proceed without splitting
# 

#Calculate DP values in R before further processing

#- locus DP (these have been updated to Nd from the Nf copy over)
## from (Nc_variant_quality_exploratory.r)
#
#    libOne  
#max (min is 0)   2152.4
#mean    15.61641   
#median  16.4
#sd  15.58382
#mean + 2*sd    46.78404
#mean - 1*sd 0.0325912


#We will aplly a mean + 2*SD filter for max DP. Given that mean - 2*SD is below 0 for all libs and mean - 1*SD is 3 or below in all but one case, let's look at distribution of DP vs GQ. The SD based filters make sense for filtering out structural variants , e.g., repetitive regions with low unique mappings (Lou et al. 2021), but if we are trying to filter out low confidence calls, we could just use quality... the only problem being that is *very* common to apply minDP filters. Depending how the DP v GQ looks, maybe just filter at DP >= 3 and call it a day. *Note that we are planning to apply and RGQ filter of >= 30 so it would be consistent to filter based on GQ*
#let's go with a minimum DP of 2 and minimum GQ of 30. Because we are haploid we have greater confidence as well... i.e., we do not need to predict/worry about heterozygotes (Lou et al. 2021)
#We will also pick up low confidence SNPs with percent NA filters (Lou et al. 2021)

#we just keep the lib one prefix here to make things easier as we move forward
vcftools --vcf INFOfilters.removed.recode.vcf --out DPGQ_filter --recode --minDP 2 --minGQ 30 --maxDP 46

#count number flagged missing
grep -v "^#" INFOfilters.removed.recode.vcf | wc -l

#488725 * 5 = 2443625
grep -v "^#" INFOfilters.removed.recode.vcf | grep -o "\s\.:" | wc -l
#203061
grep -v "^#" DPGQ_filter.recode.vcf | grep -o "\s\.\/\.:" | wc -l
#268801
#268801-203061 = 65740
# 65740/2443625 = 0.0.02690265
# 268801/2443625 = 0.1100009


#missing data filters after marking sites missing based on depth
#iterate on the cutoff vars (lmiss == locus; imiss == individual)

# the LYLF paper applies: missing inds per genotype (geno) < 50%, missing data 
# individuals (imiss) < 90%; geno < 40%, imiss < 70%; geno < 30%, imiss < 50%; 
# geno < 5%, imiss < 25%


#get percent missing per site and indv to lok at distribution to identify
# levels for filtering
# # need to sub back in . for ./.
# Be careful because vcftools rearranges the FORMAT order, BUT GT is still the 
# first value. If we want to do other manipulations we will need to take the
# FORMAT order into account
sed 's@\./\.@\.@g' DPGQ_filter.recode.vcf > DPGQ_filter.hap.vcf

vcftools --vcf DPGQ_filter.hap.vcf --missing-indv
vcftools --vcf DPGQ_filter.hap.vcf --missing-site

mv out.imiss DPGQ_filter.imiss
mv out.lmiss DPGQ_filter.lmiss
#that fixed it. The "N_data" field now has the correct number of obs, i.e., 24
# the number of samples. The ploidy issue apparently did not affect the imiss


###########################################
#First pass -- all individuals are below 15% (actually below 12%)
# apply a single hard filter to the SNPs to 25% (so max of 1 of 5 inds missing)

#somehow a comma jumped back in here. It's in the SF INFO field and looks like
#must have been inserted by vcftools (eyeroll)
perl -pe 's/(".+?[^\\]")/($ret = $1) =~ (s#,##g); $ret/ge' DPGQ_filter.hap.vcf > DPGQ_filter.hap.noComma.vcf

lmiss=0.25

awk -v var="$lmiss" '$6 > var' DPGQ_filter.lmiss | cut -f1,2 > NA.${lmiss}.sites
tail -n +2 NA.${lmiss}.sites > temp.sites && mv temp.sites NA.${lmiss}.sites
wc -l NA.${lmiss}.sites
wc -l DPGQ_filter.lmiss
#488725-70947 = 417778
#70947/488725 = 0.14516
#removing loci with >25% missing sites (70947 loci; 488725 total; 14.5%)
vcftools --vcf DPGQ_filter.hap.noComma.vcf --exclude-positions NA.${lmiss}.sites --recode --out NA_filter.loc_gt_$lmiss
# kept 417778 out of a possible 488725 Sites

#Recalculate missingness
vcftools --vcf NA_filter.loc_gt_$lmiss.recode.vcf --missing-indv
vcftools --vcf NA_filter.loc_gt_$lmiss.recode.vcf --missing-site

mv out.imiss NA_filter.loc_gt_$lmiss.imiss
mv out.lmiss NA_filter.loc_gt_$lmiss.lmiss

##########################
#cleanup and archive
#
#Final table (pass to filtering for polyalleles, MAC, LD as necessary depending on analysis)
cd ~/repo/neonectria_SNP/data/Nc
cp NA_filter.loc_gt_0.25.recode.vcf FINAL_snp.vcf
bcftools view -Ob -o FINAL_snp.bcf FINAL_snp.vcf &

#Final list of sites/indvs to filter from invairant sites tables
cat NA.*.indv > all_to_filter.indvx
cat NA.*.sites > all_missing_filter.sites
#move intermediate files and zip
mkdir filtering_intermediates
mv NA* filtering_intermediates

mv INFOfilters* filtering_intermediates

tar -czvf filtering_intermediates.tar.gz filtering_intermediates &
tar -cvf filtering_intermediates_gz.tar filtering_intermediates_gz &
rm -r filtering_intermediates
rm -r filtering_intermediates_gz

######################################

######################################

######################################
######################################
######################################
#Invariant sites filtration
#
# 1. GATK filters; remove all non-PASS
# 2. split by lib and apply genotype-level DP (based on variant site means per lib), GQ, and RGQ filters
# 3. remove lists of indv and sites from variant filtration above
# 4. apply max 25% missing filter again (will only hit ref sites; make sure to retain alist and count number filtered)
#
# Note: we may want to correct/account for the number of variant/invariant sites removed by missingness filters
# when comparing div metrics between species. Possible that there is some effect of assembly quality on the number of sites called.
# Maybe chi-square test
#
# The invariants site vcf is on the server...

# pull back out the lists of sites and indv for filtering
#tar -zxvf filtering_intermediates.tar.gz filtering_intermediates/NA*sites &
#tar -zxvf filtering_intermediates.tar.gz filtering_intermediates/NA*indv &
#also get list of sites filtered by INFO filters 
tar -zxvf filtering_intermediates.tar.gz filtering_intermediates/INFOfilters.noComma.vcf
grep -v "^#" filtering_intermediates/INFOfilters.noComma.vcf | grep -v "PASS" | cut -f 1,2 > INFO_filters.sites

wc -l INFO_filters.sites
#87050

cat all_missing_filter.sites INFO_filters.sites > all_to_filter.sites
wc -l all_to_filter.sites
#157997
#157997 - 87050 = 70947 # the difference is the amount filtered based on missingness


#remove sites and indvs (bcftools is much faster for this than vcftools) prepending with ^ means negative filter
##there are no indv to filter so we removed `-S ^all_to_filter.indv` from the command 
bcftools view -T ^all_to_filter.sites -Oz -o SNPs_INFO_and_missing_filtered.invariant_sites.vcf.gz out.invariant_sites.vcf.gz 

gunzip -c SNPs_INFO_and_missing_filtered.invariant_sites.vcf.gz | grep -v "^#" | wc -l
#42584909
gunzip -c out.invariant_sites.vcf.gz | grep -v "^#" | wc -l
# 42742906 - 42584909 = 157997
gunzip -c SNPs_INFO_and_missing_filtered.invariant_sites.vcf.gz | less -S
#32 samples confirmed

bcftools view SNPs_INFO_and_missing_filtered.invariant_sites.vcf.gz -Ob -o SNPs_INFO_and_missing_filtered.invariant_sites.bcf 

# need to filter lib ID files for the samples we have already removed
grep -vf all_to_filter.indv lib_one_ids.txt > lib_one_ids.filtered.txt
grep -vf all_to_filter.indv lib_two_ids.txt > lib_two_ids.filtered.txt
grep -vf all_to_filter.indv lib_three_ids.txt > lib_three_ids.filtered.txt
grep -vf all_to_filter.indv lib_four_ids.txt > lib_four_ids.filtered.txt

#split by library before applying depth based filters
bcftools view -S filtering_lists/lib_one_ids.filtered.txt -Ob -o lib_one.invariant.bcf SNPs_INFO_and_missing_filtered.invariant_sites.bcf &
bcftools view -S filtering_lists/lib_two_ids.filtered.txt -Ob -o lib_two.invariant.bcf SNPs_INFO_and_missing_filtered.invariant_sites.bcf &
bcftools view -S filtering_lists/lib_three_ids.filtered.txt -Ob -o lib_three.invariant.bcf SNPs_INFO_and_missing_filtered.invariant_sites.bcf & 
bcftools view -S filtering_lists/lib_four_ids.filtered.txt -Ob -o lib_four.invariant.bcf SNPs_INFO_and_missing_filtered.invariant_sites.bcf &

#filtering
bcftools +setGT lib_one.invariant.bcf -Ob -o lib_one.invariant.DP-GQ_filter.bcf -- -t q -n . -i 'FMT/DP<2 | FMT/DP>24 | GQ<30 | RGQ<30' &
bcftools +setGT lib_two.invariant.bcf -Ob -o lib_two.invariant.DP-GQ_filter.bcf -- -t q -n . -i 'FMT/DP<2 | FMT/DP>63 | GQ<30 | RGQ<30' &
bcftools +setGT lib_three.invariant.bcf -Ob -o lib_three.invariant.DP-GQ_filter.bcf -- -t q -n . -i 'FMT/DP<2 | FMT/DP>103 | GQ<30 | RGQ<30' &
bcftools +setGT lib_four.invariant.bcf -Ob -o lib_four.invariant.DP-GQ_filter.bcf -- -t q -n . -i 'FMT/DP<2 | FMT/DP>78 | GQ<30 | RGQ<30' &

#merging
for i in one two three four
do(
    bcftools index lib_${i}.invariant.DP-GQ_filter.bcf &
)
done

bcftools merge -Ob -o all_libs.invariant.DP-GQ_filtered.bcf lib_one.invariant.DP-GQ_filter.bcf lib_two.invariant.DP-GQ_filter.bcf lib_three.invariant.DP-GQ_filter.bcf lib_four.invariant.DP-GQ_filter.bcf

#filter sites on missingness > 25%
bcftools view -e 'F_MISSING>0.25' -Ob -o FINAL_invariant.bcf all_libs.invariant.DP-GQ_filtered.bcf

#check number of sites removed 
bcftools view -Ov FINAL_invariant.bcf | grep -v "^#" | wc -l
#38538664
bcftools view -Ov all_libs.invariant.DP-GQ_filtered.bcf | grep -v "^#" | wc -l
# should be 43924891 methinks; it is
# 43924891 - 38538664 = 5386227
#
bcftools view -Ov FINAL_invariant.bcf | less -S
#count hits to mt genome
bcftools view FINAL_invariant.bcf | grep -E '^LDPL01000025.1|^LDPL01000147.1|^LDPL01000119.1|^LDPL01000155.1|^LDPL01000157.1|^LDPL01000133.1|^LDPL01000153.1|^LDPL01000084.1|^LDPL01000109.1'  | wc -l
#3510 remove these
bcftools index FINAL_invariant.bcf
bcftools view FINAL_invariant.bcf | grep -v -E '^LDPL01000025.1|^LDPL01000147.1|^LDPL01000119.1|^LDPL01000155.1|^LDPL01000157.1|^LDPL01000133.1|^LDPL01000153.1|^LDPL01000084.1|^LDPL01000109.1' | bcftools view -Ob -o FINAL_invariant.nuclear.bcf
bcftools view FINAL_invariant.nuclear.bcf | grep -E '^LDPL01000025.1|^LDPL01000147.1|^LDPL01000119.1|^LDPL01000155.1|^LDPL01000157.1|^LDPL01000133.1|^LDPL01000153.1|^LDPL01000084.1|^LDPL01000109.1'  | wc -l


#Once the filtering is complete we can remove the filtering_intermediates folder
rm -R filtering_intermediates
#cleanup intermediates
mkdir filtering_intermediates_invariant
mv lib*bcf* filtering_intermediates_invariant
mv all_libs*filtered.bcf filtering_intermediates_invariant
mv SNPs_INFO* filtering_intermediates_invariant
tar -cvf filtering_intermediates_invariant.tar filtering_intermediates_invariant
rm -R filtering_intermediates_invariant

mkdir unfiltered_vcfs
mv out* unfiltered_vcfs
bgzip unfiltered_vcfs/out.vcf
tar -cvf unfiltered_vcfs.tar unfiltered_vcfs
rm -R unfiltered_vcfs

mkdir filtering_lists
mv lib*txt filtering_lists
mv *sites filtering_lists
mv *indv filtering_lists

#clean up final tables
mkdir final_tables 
mv *.bcf final_tables


#######################
#######################

#final filters for analyses
cd final_tables

#count mt hits
grep -E '^LDPL01000025.1|^LDPL01000147.1|^LDPL01000119.1|^LDPL01000155.1|^LDPL01000157.1|^LDPL01000133.1|^LDPL01000153.1|^LDPL01000084.1|^LDPL01000109.1' FINAL_snp.vcf | wc -l
#660 remove these
grep -v -E '^LDPL01000025.1|^LDPL01000147.1|^LDPL01000119.1|^LDPL01000155.1|^LDPL01000157.1|^LDPL01000133.1|^LDPL01000153.1|^LDPL01000084.1|^LDPL01000109.1' FINAL_snp.vcf > FINAL_snp.nuclear.vcf
grep -E '^LDPL01000025.1|^LDPL01000147.1|^LDPL01000119.1|^LDPL01000155.1|^LDPL01000157.1|^LDPL01000133.1|^LDPL01000153.1|^LDPL01000084.1|^LDPL01000109.1' FINAL_snp.nuclear.vcf | wc -l
bcftools view -Ob -o FINAL_snp.nuclear.bcf FINAL_snp.nuclear.vcf

#there are some samples of different individuals from the same tree or the same indv sequenced twice;
# NG163, NG20 (MI1 1.1.1); NG27, NG144 (MI1 2);
#look at which ones have better completeness
cd ~/repo/neonectria_SNP/data/Nc/final_tables/
vcftools --vcf FINAL_snp.nuclear.vcf --missing-indv
grep -E 'NG163|NG20|NG27|NG144' out.imiss

#remove dups from vcfs as well
printf 'NG20\nNG27' > ../filtering_lists/dups.txt
mkdir rm_dups

less ../filtering_lists/dups.txt

bcftools view -S ^../filtering_lists/dups.txt -Oz -o rm_dups/FINAL_invariant.nuclear.vcf.gz FINAL_invariant.nuclear.bcf &
bcftools view -S ^../filtering_lists/dups.txt -Oz -o rm_dups/FINAL_snp.IBD_analyses.vcf.gz FINAL_snp.nuclear.bcf &
bcftools view rm_dups/FINAL_invariant.nuclear.vcf.gz | grep -v '^#' | wc -l
#38535154
# this is kind of low... Nf is only about 1M off of the full genome size. Will want to look at calcs with the correction based on full genome size and the number of invariant sites as well
bcftools view rm_dups/FINAL_snp.IBD_analyses.vcf.gz | grep -v '^#' | wc -l
#1599656

#rm singletons
cd rm_dups
vcftools --gzvcf FINAL_snp.IBD_analyses.vcf.gz --mac 2 --recode --out FINAL_snp.mac_ge2
#kept 30 out of 30 Individuals
# 488739 out of a possible 1599656 Sites
bgzip FINAL_snp.mac_ge2.recode.vcf
#
#LD filter
bcftools +prune -m 0.5 -w 10000 -Oz -o FINAL_snp.mac_ge2.LD.pca_analyses.vcf.gz FINAL_snp.mac_ge2.recode.vcf.gz &

#
#bialleles & LD
bcftools view -m2 -M2 -v snps FINAL_snp.mac_ge2.recode.vcf.gz | bcftools +prune -m 0.5 -w 10000 -Oz -o FINAL_snp.mac_ge2.biallele.LD.structure_analyses.vcf.gz &

#
#bialleles
bcftools view -m2 -M2 -v snps FINAL_snp.mac_ge2.recode.vcf.gz -Oz -o FINAL_snp.mac_ge2.biallele.gwas_analyses.vcf.gz &

########################
# After final filtration we filter the metadata to only the samples retained
#
# metadata_collate/filter_to_retained_samples.R
