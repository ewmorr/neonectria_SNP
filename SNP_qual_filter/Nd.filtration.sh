
#run locally
conda activate bcftools

cd ~/repo/neonectria_SNP/data/Nd/

QDparam=2
MQparam=40
FSparam=60

gatk VariantFiltration -R ref.fasta -V out.vcf -O INFOfilters.vcf \
    --filter-name "MQ40" -filter "MQ < 40.0" \
    --filter-name "QD2" -filter "QD < 2.0" \
    --filter-name "FS60" -filter "FS > 60.0"
    

#remove commas from vcf INFO field descriptions before proceeding with vcftools
# https://unix.stackexchange.com/questions/48672/only-remove-commas-embedded-within-quotes-in-a-comma-delimited-file

perl -pe 's/(".+?[^\\]")/($ret = $1) =~ (s#,##g); $ret/ge' INFOfilters.vcf > INFOfilters.noComma.vcf


# For the diversity analyses we will keep all mac and all allele counts. We have other qual fiulters to depend on to filter out bad calls, and these types of calls are obviously important to calculation of pi and Tajima's D, for example. Also see Lou et al 2021 which points out these should not be removed for div analyses

vcftools --vcf INFOfilters.noComma.vcf --out INFOfilters.removed --recode --remove-filtered-all
# After filtering, kept 1785397 out of a possible 1983461 Sites
# 1983461 - 1785397 = 198064

#There are no extreme %NA outliers so we just proceed with filtering as normal (as opposed to Nf)

####################
#DP filters
# we are prforming filtering on the INFO filed filtered dataset, but NOT having removed biallelic and MAC filters.
# This dataset will be used for diversity metrics. These are presumably good quality calls based on the other filters
#For any GWAS or GxE and STRUCTURE we will want to filter for biallele and MAC as the methods either ignore or
# don't handle unfiltered
# can perform biallelic and MAC  filters on the table after performing the DP and NA filters

#note we tested for library effects on singleton alleles per individual and 
# found no significant differences (see Nf_MAC_per_indv.r and 
# Nf/MAC_singletons_per_indv_by_library.pdf) so even if these are errorneous 
# they are randomly distributed and not a technical artefact.

# we first need to split the dataset by library
#library IDs are generated with write_library_sample_ids.r

vcftools --vcf INFOfilters.removed.recode.vcf --out lib_one --recode --keep lib_one_ids.txt
vcftools --vcf INFOfilters.removed.recode.vcf --out lib_two --recode --keep lib_two_ids.txt
vcftools --vcf INFOfilters.removed.recode.vcf --out lib_three --recode --keep lib_three_ids.txt
vcftools --vcf INFOfilters.removed.recode.vcf --out lib_four --recode --keep lib_four_ids.txt

#Calculate DP values in R before further processing

#- locus DP across libraries (these have been updated to Nd from the Nf copy over)
## from (Nd_variant_quality_exploratory.r)
#
#    libOne  libTwo  libThree    libFour
#max (min is 0)   943.625 2978.222  6253   3964.214
#mean    9.332779    26.79611   66.00215  39.38802
#median  9.625    29.66667   75  43.57143
#sd  7.480772 16.13664    38.98371    29.10379
#mean + 2*sd    24.29432    59.06939    143.9696   97.59561
#mean - 1*sd 1.852007   10.65946 27.01844   10.28423


#We will aplly a mean + 2*SD filter for max DP. Given that mean - 2*SD is below 0 for all libs and mean - 1*SD is 3 or below in all but one case, let's look at distribution of DP vs GQ. The SD based filters make sense for filtering out structural variants , e.g., repetitive regions with low unique mappings (Lou et al. 2021), but if we are trying to filter out low confidence calls, we could just use quality... the only problem being that is *very* common to apply minDP filters. Depending how the DP v GQ looks, maybe just filter at DP >= 3 and call it a day. *Note that we are planning to apply and RGQ filter of >= 30 so it would be consistent to filter based on GQ*
#let's go with a minimum DP of 2 and minimum GQ of 30. Because we are haploid we have greater confidence as well... i.e., we do not need to predict/worry about heterozygotes (Lou et al. 2021)
#We will also pick up low confidence SNPs with percent NA filters (Lou et al. 2021)


vcftools --vcf lib_one.recode.vcf --out lib_one.DPGQ_filter --recode --minDP 2 --minGQ 30 --maxDP 24
vcftools --vcf lib_two.recode.vcf --out lib_two.DPGQ_filter --recode --minDP 2 --minGQ 30 --maxDP 59
vcftools --vcf lib_three.recode.vcf --out lib_three.DPGQ_filter --recode --minDP 2 --minGQ 30 --maxDP 143
vcftools --vcf lib_four.recode.vcf --out lib_four.DPGQ_filter --recode --minDP 2 --minGQ 30 --maxDP 97

#count number flagged missing
grep -v "^#" lib_one.recode.vcf | wc -l
wc -l lib_one_ids.txt
#1785397 * 8 = 14283176
grep -v "^#" lib_one.recode.vcf | grep -o "\s\.:" | wc -l
#1083904
grep -v "^#" lib_one.DPGQ_filter.recode.vcf | grep -o "\s\.\/\.:" | wc -l
#1451645
#1451645-1083904 = 367741
# 367741/14283176 = 0.02574644
# 1451645/14283176 = 0.1016332

#grep -v "^#" lib_two.recode.vcf | wc -l #all have the same number of SNPs
wc -l lib_two_ids.txt
# 1322180*9 = 11899620
grep -v "^#" lib_two.recode.vcf | grep -o "\s\.:" | wc -l
#1048464
grep -v "^#" lib_two.DPGQ_filter.recode.vcf | grep -o "\s\.\/\.:" | wc -l
#2439573
#2439573-1048464=1391109
#1048464/11899620 = 0.08810903
# 2439573/11899620 = 0.2050127

#######################
#These figures need to be updated
#for lib3 and 4
wc -l lib_three_ids.txt
#6*1322180 = 7933080
grep -v "^#" lib_three.recode.vcf | grep -o "\s\.:" | wc -l
#546987
grep -v "^#" lib_three.DPGQ_filter.recode.vcf | grep -o "\s\.\/\.:" | wc -l
#790064 - 546987 = 243077
# 243077/7933080 = 0.03064094
# 790064/7933080 = 0.09959108

wc -l lib_four_ids.txt
# 11*1322180 = 14543980
grep -v "^#" lib_four.recode.vcf | grep -o "\s\.:" | wc -l
#1014567
grep -v "^#" lib_four.DPGQ_filter.recode.vcf | grep -o "\s\.\/\.:" | wc -l
#1323941 - 1014567 = 309374
# 309374/14543980 = 0.02127162
# 1323941/14543980 = 0.09103017
###################################

#missing data filters after marking sites missing based on depth
#iterate on the cutoff vars (lmiss == locus; imiss == individual)

# the LYLF paper applies: missing inds per genotype (geno) < 50%, missing data 
# individuals (imiss) < 90%; geno < 40%, imiss < 70%; geno < 30%, imiss < 50%; 
# geno < 5%, imiss < 25%

#first merge the vcfs from the four libs.
# need to bgzip and tabix first

for i in one two three four
do(
    bgzip lib_${i}.DPGQ_filter.recode.vcf
    tabix lib_${i}.DPGQ_filter.recode.vcf.gz
)
done

vcf-merge lib_one.DPGQ_filter.recode.vcf.gz \
    lib_two.DPGQ_filter.recode.vcf.gz \
    lib_three.DPGQ_filter.recode.vcf.gz \
    lib_four.DPGQ_filter.recode.vcf.gz > all_libs.DPGQ_filter.recode.vcf
    
#get percent missing per site and indv to lok at distribution to identify
# levels for filtering
vcftools --vcf all_libs.DPGQ_filter.recode.vcf --missing-indv
vcftools --vcf all_libs.DPGQ_filter.recode.vcf --missing-site
# interesting... looking at the result of above, the lmiss values are incorrect
# i.e., it looks like vcftools is counting ./. (that it inserted...) as two 
# observations, whereas good calls are 1 obs (they should all be one), Try 
# subbing . back in for ./. and recalculate

# Be careful because vcftools rearranges the FORMAT order, BUT GT is still the 
# first value. If we want to do other manipulations we will need to take the
# FORMAT order into account
sed 's@\./\.@\.@g' all_libs.DPGQ_filter.recode.vcf > all_libs.DPGQ_filter.hap.vcf

vcftools --vcf all_libs.DPGQ_filter.hap.vcf --missing-indv
vcftools --vcf all_libs.DPGQ_filter.hap.vcf --missing-site

mv out.imiss all_libs.DPGQ_filter.imiss
mv out.lmiss all_libs.DPGQ_filter.lmiss
#that fixed it. The "N_data" field now has the correct number of obs, i.e., 24
# the number of samples. The ploidy issue apparently did not effect the imiss


###########################################
#First pass -- there is one individual that is above 60% missingness after 
# applying the DP and GQ filters, whereas all others are below 15%
# let's remove this (so everything now below 15% like Nf) and then apply a 
# single hard filter to the SNPs to 25%


imiss=0.15

awk -v var="$imiss" '$5 > var' all_libs.DPGQ_filter.imiss | cut -f1 > NA.${imiss}.indv
tail -n +2 NA.${imiss}.indv > temp.sites && mv temp.sites NA.${imiss}.indv


#somehow a comma jumped back in here. It's in the SF INFO field and looks like
#must have been inserted by vcftools (eyeroll)
perl -pe 's/(".+?[^\\]")/($ret = $1) =~ (s#,##g); $ret/ge' all_libs.DPGQ_filter.hap.vcf > all_libs.DPGQ_filter.hap.noComma.vcf

vcftools --vcf all_libs.DPGQ_filter.hap.noComma.vcf --remove NA.${imiss}.indv --recode --out NA_filter.ind_gt_${imiss}

#Recalculate missingness
vcftools --vcf NA_filter.ind_gt_${imiss}.recode.vcf --missing-indv
vcftools --vcf NA_filter.ind_gt_${imiss}.recode.vcf --missing-site

mv out.imiss NA_filter.ind_gt_${imiss}.imiss
mv out.lmiss NA_filter.ind_gt_${imiss}.lmiss


# 
# ##2nd pass

#
## ##3rd pass

lmiss=0.25

awk -v var="$lmiss" '$6 > var' NA_filter.ind_gt_${imiss}.lmiss | cut -f1,2 > NA.${lmiss}.sites
tail -n +2 NA.${lmiss}.sites > temp.sites && mv temp.sites NA.${lmiss}.sites
#removing loci with >25% missing sites (185081 loci; 1785397 total; 10.4%)
vcftools --vcf NA_filter.ind_gt_${imiss}.recode.vcf --exclude-positions NA.${lmiss}.sites --recode --out NA_filter.loc_gt_$lmiss.ind_gt_${imiss}

#Recalculate missingness
vcftools --vcf NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.recode.vcf --missing-indv
vcftools --vcf NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.recode.vcf --missing-site

mv out.imiss NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.imiss
mv out.lmiss NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.lmiss

##########################
#cleanup and archive
#
#Final table (pass to filtering for polyalleles, MAC, LD as necessary depending on analysis)
cd ~/repo/neonectria_SNP/data/Nd
cp NA_filter.loc_gt_0.25.ind_gt_0.15.recode.vcf FINAL_snp.vcf
bcftools view -Ob -o FINAL_snp.bcf FINAL_snp.vcf &

#Final list of sites/indvs to filter from invairant sites tables
cat NA.*.indv > all_to_filter.indv
cat NA.*.sites > all_missing_filter.sites
#move intermediate files and zip
mkdir filtering_intermediates
mv NA* filtering_intermediates
mv lib*vcf filtering_intermediates
ls lib*
mv lib*log filtering_intermediates

mv all_libs* filtering_intermediates
mv INFOfilters* filtering_intermediates

mkdir filtering_intermediates_gz
mv lib*gz* filtering_intermediates_gz

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
#198064

cat all_missing_filter.sites INFO_filters.sites > all_to_filter.sites
wc -l all_to_filter.sites
#383145
#383145 - 198064 = 185081 # the difference is the amount filtered based on missingness


#remove sites and indvs (bcftools is much faster for this than vcftools) prepending with ^ means negative filter
bcftools view -T ^all_to_filter.sites -S ^all_to_filter.indv out.invariant_sites.vcf.gz | bgzip > SNPs_INFO_and_missing_filtered.invariant_sites.vcf.gz

gunzip -c SNPs_INFO_and_missing_filtered.invariant_sites.vcf.gz | grep -v "^#" | wc -l
#43924891
gunzip -c out.invariant_sites.vcf.gz | grep -v "^#" | wc -l
# 44308036 - 43924891 = 383145
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
cd ~/repo/neonectria_SNP/data/Nd/final_tables/
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
