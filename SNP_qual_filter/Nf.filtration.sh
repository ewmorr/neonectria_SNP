
#run locally
conda activate bcftools

cd ~/repo/neonectria_SNP/data/Nf/

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

#remove >2 alleles, mac<2, all non-PASS

vcftools --vcf INFOfilters.noComma.vcf --out INFOfilters_mac2_bialllic.vcf --max-alleles 2 --mac 2 --remove-filtered-all
# After filtering, kept 620910 out of a possible 1403293 Sites

grep -v "^#" INFOfilters.noComma.vcf | wc -l
#1403293
grep -v "^#" INFOfilters.noComma.vcf | grep -v "PASS" | wc -l
#81113

# 1403293 - 659221 = 744072 total removed
# 744072 - 81113 = 662959 removed as multi alleles or mac<2 (need to check this there may be some overlap between the INFO fields and the others

vcftools --vcf INFOfilters.noComma.vcf --out INFOfilters_bialllic.vcf --max-alleles 2
# After filtering, kept 1330626 out of a possible 1403293 Sites
#1403293 - 1330626 = 72667
vcftools --vcf INFOfilters.noComma.vcf --out INFOfilters_mac2.vcf --mac 2
#After filtering, kept 686540 out of a possible 1403293 Sites
# 1403293 - 686540 = 716753

# For the diversity analyses we will keep all mac and all allele counts. We have other qual fiulters to depend on to filter out bad calls, and these types of calls are obviously important to calculation of pi and Tajima's D, for example. Also see Lou et al 2021 which points out these should not be removed for div analyses

vcftools --vcf INFOfilters.noComma.vcf --out INFOfilters.removed.vcf --remove-filtered-all
# After filtering, kept 1322180 out of a possible 1403293 Sites
# 1403293 - 1322180 = 81113

#We also filter NG130 because it has very high NA percent (~60%) relative to all  the other indvs (<20%) see individual.NA_DP.pdf

vcftools --vcf INFOfilters.noComma.vcf --out INFOfilters.removed --recode --remove-filtered-all --remove-indv "NG130"

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

#vcftools calcs aon MAC==1 for comnparison to vcfR
vcftools --vcf INFOfilters.removed.recode.vcf --out temp --mac 2
#1322180 - 646496 = 675684
vcftools --vcf INFOfilters.removed.recode.vcf --out temp --mac 1
#1322180 - 979522 = 342658
#675684 - 342658 = 333026

# we first need to split the dataset by library
#library IDs are generated with write_library_sample_ids.r

vcftools --vcf INFOfilters.removed.recode.vcf --out lib_one --recode --keep lib_one_ids.txt
vcftools --vcf INFOfilters.removed.recode.vcf --out lib_two --recode --keep lib_two_ids.txt
vcftools --vcf INFOfilters.removed.recode.vcf --out lib_three --recode --keep lib_three_ids.txt
vcftools --vcf INFOfilters.removed.recode.vcf --out lib_four --recode --keep lib_four_ids.txt

#- locus DP across libraries (Nf_variant_quality_exploratory.r)
#    libOne  libTwo  libThree    libFour
#max (min is 0)   840.7727 1840.81  2335.5   2264.636
#mean    8.1    22.8   36.5  29.1
#median  6.5    21.9   36.3  28.8
#sd  8.2 20.3    33.5    24.9
#mean + 2*sd    24.5    63.4    103.5   78.9
#mean - 1*sd -0.1   2.5 3   4.2

#filter based on (0.25 or 3x of mean)
#minDP = ceiling(mean/4)
#maxDP = floor(mean * 3)
#   libOne  libTwo  libThree    libFour
#minDP 3     6   10 8
#max    24  68  109 87
#These seem pretty arbitrary, and hopefully we can get beyond just multiplying/dividing by a random factor....

#We will apply a mean + 2*SD filter for max DP. Given that mean - 2*SD is below 0 for all libs and mean - 1*SD is 3 or below in all but one case, let's look at distribution of DP vs GQ. The SD based filters make sense for filtering out structural variants , e.g., repetitive regions with low unique mappings (Lou et al. 2021), but if we are trying to filter out low confidence calls, we could just use quality... the only problem being that is *very* common to apply minDP filters. Depending how the DP v GQ looks, maybe just filter at DP >= 3 and call it a day. *Note that we are planning to apply and RGQ filter of >= 30 so it would be consistent to filter based on GQ*
#let's go with a minimum DP of 2 and minimum GQ of 30. Because we are haploid we have greater confidence as well... i.e., we do not need to predict/worry about heterozygotes (Lou et al. 2021)
#We will also pick up low confidence SNPs with percent NA filters (Lou et al. 2021)

#filter per lib
vcftools --vcf lib_one.recode.vcf --out lib_one.DPGQ_filter --recode --minDP 2 --minGQ 30 --maxDP 24
vcftools --vcf lib_two.recode.vcf --out lib_two.DPGQ_filter --recode --minDP 2 --minGQ 30 --maxDP 63
vcftools --vcf lib_three.recode.vcf --out lib_three.DPGQ_filter --recode --minDP 2 --minGQ 30 --maxDP 103
vcftools --vcf lib_four.recode.vcf --out lib_four.DPGQ_filter --recode --minDP 2 --minGQ 30 --maxDP 78

#count number flagged missing
grep -v "^#" lib_one.recode.vcf | wc -l
wc -l lib_one_ids.txt
#1322180 * 66 = 87263880
grep -v "^#" lib_one.recode.vcf | grep -o "\s\.:" | wc -l
#7081788
grep -v "^#" lib_one.DPGQ_filter.recode.vcf | grep -o "\s\.:" | wc -l
#0
#now that vcftools has modified the GTs directly NA is recorded as './.' instead of '.'

grep -v "^#" lib_one.DPGQ_filter.recode.vcf | grep -o "\s\.\/\.:" | wc -l
#11261826
#11261826-7081788 = 4180038
# 4180038/87263880 = 0.04790112
# 11261826/87263880 = 0.1290548

#grep -v "^#" lib_two.recode.vcf | wc -l #all have the same number of SNPs
wc -l lib_two_ids.txt
# 1322180*42 = 55531560
grep -v "^#" lib_two.recode.vcf | grep -o "\s\.:" | wc -l
#3796076
grep -v "^#" lib_two.DPGQ_filter.recode.vcf | grep -o "\s\.\/\.:" | wc -l
#6464428
#6464428-3796076=2668352
#2668352/55531560 = 0.04805109
# 6464428/55531560 = 0.11641

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

#locus missingness looks approximately half normal distribution centered on zero
#2*SD = 0.4288
#95th percentile = 0.524194
#
#imiss approximately normal
#mean = 0.1210 
#mean+2*SD = 0.2424
#95th percentile = 0.2208
#
# see Nf.lmiss_imiss_dist.r and Nf.lmiss_imiss_dist.pdf
# red lines are 95th percentile and blue are mean+2*SD
# 
# We will filter top 5% of loci and reassess
#
##The following is run iteratively with Nf.lmiss_imiss_dist.pretty_plot.r
# ##1st pass

lmiss=0.524194

awk -v var="$lmiss" '$6 > var' all_libs.DPGQ_filter.lmiss | cut -f1,2 > NA.${lmiss}.sites
wc -l NA.${lmiss}.sites
wc -l all_libs.DPGQ_filter.lmiss


#somehow a comma jumped back in here. It's in the SF INFO field and looks like
#must have been inserted by vcftools (eyeroll)
perl -pe 's/(".+?[^\\]")/($ret = $1) =~ (s#,##g); $ret/ge' all_libs.DPGQ_filter.hap.vcf > all_libs.DPGQ_filter.hap.noComma.vcf

#no header line is expected in the input
tail -n +2 NA.${lmiss}.sites > temp.sites && mv temp.sites NA.${lmiss}.sites

#removing top 5% missing loci
vcftools --vcf all_libs.DPGQ_filter.hap.noComma.vcf --exclude-positions NA.${lmiss}.sites --recode --out NA_filter.loc_gt_${lmiss}

#Recalculate missingness
vcftools --vcf NA_filter.loc_gt_${lmiss}.recode.vcf --missing-indv
vcftools --vcf NA_filter.loc_gt_${lmiss}.recode.vcf --missing-site

mv out.imiss NA_filter.loc_gt_${lmiss}.imiss
mv out.lmiss NA_filter.loc_gt_${lmiss}.lmiss

#lmiss
#2*SD = 0.300844
#95th = 0.36290300 

#imiss -- there is one outlier well above 50% rest are well below 40%
#mean = 0.09133246
#mean + 2*SD = 0.2166592
#95th = 0.19588090
##99th = 0.30502884
##
#
# After removing the top 5% of locus NA there is one clear outlier in indv
# We could use 40% cutoff to remove this indv
# but let's stick with 99th percentile to be consistent
# and iterate imiss 99th lmis 95th
# 
# ##2nd pass
imiss=0.30502884

awk -v var="$imiss" '$5 > var' NA_filter.loc_gt_${lmiss}.imiss | cut -f1 > NA.${imiss}.indv
tail -n +2 NA.${imiss}.indv > temp.sites && mv temp.sites NA.${imiss}.indv
vcftools --vcf NA_filter.loc_gt_${lmiss}.recode.vcf --remove NA.${imiss}.indv --recode --out NA_filter.loc_gt_$lmiss.ind_gt_${imiss}

#Recalculate missingness
vcftools --vcf NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.recode.vcf --missing-indv
vcftools --vcf NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.recode.vcf --missing-site

mv out.imiss NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.imiss
mv out.lmiss NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.lmiss

#lmiss
#2*SD = 0.2941666
#95th = 0.3606560

#imiss 
#mean = 0.0696765
#mean + 2*SD = 0.1601633
#95th = 0.15338275
##99th = 0.24944936
##
#
## ##3rd pass

lmiss_old=$lmiss
lmiss=0.3606560
awk -v var="$lmiss" '$6 > var' NA_filter.loc_gt_$lmiss_old.ind_gt_${imiss}.lmiss | cut -f1,2 > NA.${lmiss}.sites
tail -n +2 NA.${lmiss}.sites > temp.sites && mv temp.sites NA.${lmiss}.sites
#removing top 5% missing loci
vcftools --vcf NA_filter.loc_gt_$lmiss_old.ind_gt_${imiss}.recode.vcf --exclude-positions NA.${lmiss}.sites --recode --out NA_filter.loc_gt_$lmiss.ind_gt_${imiss}

#Recalculate missingness
vcftools --vcf NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.recode.vcf --missing-indv
vcftools --vcf NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.recode.vcf --missing-site

mv out.imiss NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.imiss
mv out.lmiss NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.lmiss

#lmiss
#2*SD = 0.2300616
#95th = 0.2786890

#imiss 
#mean = 0.06837117
#mean + 2*SD = 0.149429
#95th = 0.13440880
##99th = 0.23454579
##
####99th removes top 2
####
#### ## 4th pass
#### 
imiss_old=$imiss
imiss=0.23454579

awk -v var="$imiss" '$5 > var' NA_filter.loc_gt_$lmiss.ind_gt_$imiss_old.imiss | cut -f1 > NA.${imiss}.indv
tail -n +2 NA.${imiss}.indv > temp.sites && mv temp.sites NA.${imiss}.indv
#removing top 1% missing indv, NG48, NG161
vcftools --vcf NA_filter.loc_gt_$lmiss.ind_gt_$imiss_old.recode.vcf --remove NA.${imiss}.indv --recode --out NA_filter.loc_gt_$lmiss.ind_gt_${imiss}
#Recalculate missingness 
vcftools --vcf NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.recode.vcf --missing-indv
vcftools --vcf NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.recode.vcf --missing-site

mv out.imiss NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.imiss
mv out.lmiss NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.lmiss

#lmiss
#2*SD = 0.2238579
#95th = 0.2750000

#imiss 
#mean = 0.06521228
#mean + 2*SD = 0.1300293
#95th = 0.11858750 
##99th = 0.18945936 

## 5th pass

lmiss_old=$lmiss
lmiss=0.2750000
awk -v var="$lmiss" '$6 > var' NA_filter.loc_gt_$lmiss_old.ind_gt_${imiss}.lmiss | cut -f1,2 > NA.${lmiss}.sites
tail -n +2 NA.${lmiss}.sites > temp.sites && mv temp.sites NA.${lmiss}.sites
#removing top 5% missing loci
vcftools --vcf NA_filter.loc_gt_$lmiss_old.ind_gt_${imiss}.recode.vcf --exclude-positions NA.${lmiss}.sites --recode --out NA_filter.loc_gt_$lmiss.ind_gt_${imiss}

#Recalculate missingness
vcftools --vcf NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.recode.vcf --missing-indv
vcftools --vcf NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.recode.vcf --missing-site

mv out.imiss NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.imiss
mv out.lmiss NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.lmiss

#lmiss
#2*SD = 0.1805771
#95th = 0.225000 

#imiss 
#mean = 0.05275892
#mean + 2*SD = 0.1153699
#95th = 0.10485945
##99th = 0.17720253
#99th removes top 2

## 6th pass
## 
imiss_old=$imiss
imiss=0.17720253

awk -v var="$imiss" '$5 > var' NA_filter.loc_gt_$lmiss.ind_gt_$imiss_old.imiss | cut -f1 > NA.${imiss}.indv
tail -n +2 NA.${imiss}.indv > temp.sites && mv temp.sites NA.${imiss}.indv
#removing top 1% missing indv, NG55, NG62
vcftools --vcf NA_filter.loc_gt_$lmiss.ind_gt_$imiss_old.recode.vcf --remove NA.${imiss}.indv --recode --out NA_filter.loc_gt_$lmiss.ind_gt_${imiss}
#Recalculate missingness
vcftools --vcf NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.recode.vcf --missing-indv
vcftools --vcf NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.recode.vcf --missing-site

mv out.imiss NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.imiss
mv out.lmiss NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.lmiss

#lmiss
#2*SD = 0.1761217
#95th = 0.220339 

#imiss 
#mean = 0.05040136
#mean + 2*SD = 0.1016853
#95th = 0.09642650
##99th = 0.12870211

## 7th pass

lmiss_old=$lmiss
lmiss=0.220339
awk -v var="$lmiss" '$6 > var' NA_filter.loc_gt_$lmiss_old.ind_gt_${imiss}.lmiss | cut -f1,2 > NA.${lmiss}.sites
tail -n +2 NA.${lmiss}.sites > temp.sites && mv temp.sites NA.${lmiss}.sites
#removing top 5% missing loci
vcftools --vcf NA_filter.loc_gt_$lmiss_old.ind_gt_${imiss}.recode.vcf --exclude-positions NA.${lmiss}.sites --recode --out NA_filter.loc_gt_$lmiss.ind_gt_${imiss}

#Recalculate missingness
vcftools --vcf NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.recode.vcf --missing-indv
vcftools --vcf NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.recode.vcf --missing-site

mv out.imiss NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.imiss
mv out.lmiss NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.lmiss

#lmiss
#2*SD = 0.1423921
#95th = 0.177966

#imiss 
#mean = 0.04063663
#mean + 2*SD = 0.08922933
#95th = 0.08485095
##99th = 0.11892679
##
## 8th pass

imiss_old=$imiss
imiss=0.11892679

awk -v var="$imiss" '$5 > var' NA_filter.loc_gt_$lmiss.ind_gt_$imiss_old.imiss | cut -f1 > NA.${imiss}.indv
tail -n +2 NA.${imiss}.indv > temp.sites && mv temp.sites NA.${imiss}.indv
#removing top 1% missing indv, NG55, NG62
vcftools --vcf NA_filter.loc_gt_$lmiss.ind_gt_$imiss_old.recode.vcf --remove NA.${imiss}.indv --recode --out NA_filter.loc_gt_$lmiss.ind_gt_${imiss}
#Recalculate missingness
vcftools --vcf NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.recode.vcf --missing-indv
vcftools --vcf NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.recode.vcf --missing-site

mv out.imiss NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.imiss
mv out.lmiss NA_filter.loc_gt_$lmiss.ind_gt_${imiss}.lmiss

#lmiss
#2*SD = 0.1388054
#95th = 0.17241400
##max = 0.224138

#imiss 
#mean = 0.03892257
#mean + 2*SD = 0.07986812
#95th = 0.07593268
##99th = 0.10130399
##
## max = 0.110865
## 

########################
########################
#Testing the effect of
# filtering the full list
# based on the final criteria
# 

echo $lmiss
echo $imiss

awk -v var="$lmiss" '$6 > var' all_libs.DPGQ_filter.lmiss | cut -f1,2 > test_hard_filter_NA.${lmiss}.sites
awk -v var="$imiss" '$5 > var' all_libs.DPGQ_filter.imiss | cut -f1 > test_hard_filter_NA.${imiss}.indv

wc -l test_hard_filter_NA.${lmiss}.sites
wc -l test_hard_filter_NA.${imiss}.indv

# the iterative filtering strategy removes 233959 sites and 8 indv (over the eight passes)
# whereas hard filtering with the same criteria would remove 252967-1 sites and 47-1 indv
# about 19K less sites

#how many indv removed if first filter by site hard cutoff and then by indv
vcftools --vcf all_libs.DPGQ_filter.hap.vcf --exclude-positions test_hard_filter_NA.${lmiss}.sites --recode --out test_hard_filter_NA.loc_gt_$lmiss
#recalculate indv missingness
vcftools --vcf test_hard_filter_NA.loc_gt_$lmiss.recode.vcf --missing-indv
mv out.imiss test_hard_filter_NA.loc_gt_$lmiss.imiss
awk -v var="$imiss" '$5 > var' test_hard_filter_NA.loc_gt_$lmiss.imiss | cut -f1 > test_hard_filter_NA.loc_gt_$lmiss.indv_gt_${imiss}.indv

wc -l test_hard_filter_NA.loc_gt_$lmiss.indv_gt_${imiss}.indv
#actually only 7 indv filtered this way, so one less...


########################
########################
#Testing the effect of
# filtering the full list
# based on final criteria
# of 0.25 sites and 0.15 indv

echo $lmiss
echo $imiss

awk -v var="$lmiss" '$6 > var' all_libs.DPGQ_filter.lmiss | cut -f1,2 > test_hard_filter_NA.${lmiss}.sites
awk -v var="$imiss" '$5 > var' all_libs.DPGQ_filter.imiss | cut -f1 > test_hard_filter_NA.${imiss}.indv

wc -l test_hard_filter_NA.${lmiss}.sites
#214570 sites
vcftools --vcf all_libs.DPGQ_filter.hap.vcf --exclude-positions test_hard_filter_NA.${lmiss}.sites --recode --out test_hard_filter_NA.loc_gt_$lmiss
# keeps 1107611

#recalculate indv missingness
vcftools --vcf test_hard_filter_NA.loc_gt_$lmiss.recode.vcf --missing-indv
mv out.imiss test_hard_filter_NA.loc_gt_$lmiss.imiss
awk -v var="$imiss" '$5 > var' test_hard_filter_NA.loc_gt_$lmiss.imiss | cut -f1 > test_hard_filter_NA.loc_gt_$lmiss.indv_gt_${imiss}.indv

wc -l test_hard_filter_NA.loc_gt_$lmiss.indv_gt_${imiss}.indv
#would remove 7 indv



########################
########################
#Testing the effect of
# applying 25% missing locus
# and 15% missing indv
# hard filter after the
# iterative percentile cutoffs
# filtering the full list
# based on the final criteria
# 

######################################
#THIS IS NOW THE FINAL RUN
lmiss=0.25
imiss=0.15

awk -v var="$lmiss" '$6 > var' NA_filter.loc_gt_0.2750000.ind_gt_0.17720253.lmiss | cut -f1,2 > NA.${lmiss}.sites
tail -n +2 NA.${lmiss}.sites > temp.sites && mv temp.sites NA.${lmiss}.sites
vcftools --vcf NA_filter.loc_gt_0.2750000.ind_gt_0.17720253.recode.vcf --exclude-positions NA.${lmiss}.sites --recode --out NA_filter.iterative_plus_hard_filter.loc_gt_${lmiss}.ind_gt_0.17720253
# keeps 1116217 sites

#recalc missingness
vcftools --vcf NA_filter.iterative_plus_hard_filter.loc_gt_${lmiss}.ind_gt_0.17720253.recode.vcf --missing-indv
vcftools --vcf NA_filter.iterative_plus_hard_filter.loc_gt_${lmiss}.ind_gt_0.17720253.recode.vcf --missing-site

mv out.imiss NA_filter.iterative_plus_hard_filter.loc_gt_${lmiss}.ind_gt_0.17720253.imiss
mv out.lmiss NA_filter.iterative_plus_hard_filter.loc_gt_${lmiss}.ind_gt_0.17720253.lmiss

imiss=0.15
awk -v var="$imiss" '$5 > var' NA_filter.iterative_plus_hard_filter.loc_gt_${lmiss}.ind_gt_0.17720253.imiss | cut -f1 > NA.${imiss}.indv
tail -n +2 NA.${imiss}.indv > temp.sites && mv temp.sites NA.${imiss}.indv
vcftools --vcf test_hard_filter_NA.loc_gt_${lmiss}.ind_gt_0.17720253.recode.vcf --remove NA.${imiss}.indv --recode --out NA_filter.iterative_plus_hard_filter.loc_gt_$lmiss.ind_gt_${imiss}

vcftools --vcf NA_filter.iterative_plus_hard_filter.loc_gt_$lmiss.ind_gt_${imiss}.recode.vcf --missing-indv
vcftools --vcf NA_filter.iterative_plus_hard_filter.loc_gt_$lmiss.ind_gt_${imiss}.recode.vcf --missing-site

mv out.imiss NA_filter.iterative_plus_hard_filter.loc_gt_${lmiss}.ind_gt_$imiss.imiss
mv out.lmiss NA_filter.iterative_plus_hard_filter.loc_gt_${lmiss}.ind_gt_$imiss.lmiss

#THIS IS NOW THE FINAL RUN
##########################
#cleanup and archive
#
#Final table (pass to filtering for polyalleles, MAC, LD as necessary depending on analysis)
cp NA_filter.iterative_plus_hard_filter.loc_gt_$lmiss.ind_gt_${imiss}.recode.vcf FINAL_snp.vcf

#Final list of sites/indvs to filter from invairant sites tables
cat NA.*.indv > indv_to_remove.indv
cat NA.*.sites > sites_to_exclude.sites
#move intermediate files and zip
mkdir filtering_intermediates
mv NA* filtering_intermediates
mv lib*vcf filtering_intermediates
ls lib*
mv lib*log filtering_intermediates
mv 4th* filtering_intermediates
mv all_libs* filtering_intermediates
mv INFOfilters* filtering_intermediates

mkdir filtering_intermediates_gz
mv lib*gz* filtering_intermediates_gz

tar -czvf filtering_intermediates.tar.gz filtering_intermediates
tar -cvf filtering_intermediates_gz.tar filtering_intermediates_gz
rm -r filtering_intermediates
rm -r filtering_intermediates_gz

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
tar -zxvf filtering_intermediates.tar.gz filtering_intermediates/NA*sites
tar -zxvf filtering_intermediates.tar.gz filtering_intermediates/NA*indv
#note we also made a file for NG130 which was filtered out early

#also get list of sites filtered by INFO filters 
tar -zxvf filtering_intermediates.tar.gz filtering_intermediates/INFOfilters.noComma.vcf
grep -v "^#" filtering_intermediates/INFOfilters.noComma.vcf | grep -v "PASS" | cut -f 1,2 > INFO_filters.sites

wc -l filtering_intermediates/INFO_filters.sites
#81113

cat filtering_intermediates/*sites > all_to_filter.sites
wc -l filtering_intermediates/all_to_filter.sites
#287076
#287076 - 81113 = 205963

#note that we created a file with only NG130 listed, which was removed 
# before the iterative missingness filtration bc of very high missingness
# prior to marking sites missing for DP or GQ
echo "NG130" > filtering_intermediates/NA.NG130.indv
cat filtering_intermediates/*indv > all_to_filter.indv
wc -l filtering_intermediates/all_to_filter.indv
#8
#

#remove sites and indvs (bcftools is much faster for this than vcftools) prepending with ^ means negative filter
bcftools view -T ^filtering_intermediates/all_to_filter.sites -S ^filtering_intermediates/all_to_filter.indv out.invariant_sites.vcf.gz | bgzip > SNPs_INFO_and_missing_filtered.invariant_sites.vcf.gz

gunzip -c SNPs_INFO_and_missing_filtered.invariant_sites.vcf.gz | grep -v "^#" | wc -l
#42335985
gunzip -c out.invariant_sites.vcf.gz | grep -v "^#" | wc -l
# 42623061 - 42335985 = 287076
gunzip -c SNPs_INFO_and_missing_filtered.invariant_sites.vcf.gz | less -S
gunzip -c SNPs_INFO_and_missing_filtered.invariant_sites.vcf.gz | grep -m 1  "^#CHROM" | grep -o "NG" | wc -l
#117 samples confirmed

bcftools view SNPs_INFO_and_missing_filtered.invariant_sites.vcf.gz -Ob -o SNPs_INFO_and_missing_filtered.invariant_sites.bcf

# need to filter lib ID files for the samples we have already removed
grep -wvf filtering_lists/all_to_filter.indv filtering_lists/lib_one_ids.txt > filtering_lists/lib_one_ids.filtered.txt
grep -wvf filtering_lists/all_to_filter.indv filtering_lists/lib_two_ids.txt > filtering_lists/lib_two_ids.filtered.txt
grep -wvf filtering_lists/all_to_filter.indv filtering_lists/lib_three_ids.txt > filtering_lists/lib_three_ids.filtered.txt
grep -wvf filtering_lists/all_to_filter.indv filtering_lists/lib_four_ids.txt > filtering_lists/lib_four_ids.filtered.txt

cat filtering_lists/lib*ids.filtered.txt | wc -l
#117 needed the -w optyion for exact match

#split by library before applying depth based filters
bcftools view -S filtering_lists/lib_one_ids.filtered.txt -Ob -o lib_one.invariant.bcf filtering_intermediates_invariant/SNPs_INFO_and_missing_filtered.invariant_sites.bcf &
bcftools view -S filtering_lists/lib_two_ids.filtered.txt -Ob -o lib_two.invariant.bcf filtering_intermediates_invariant/SNPs_INFO_and_missing_filtered.invariant_sites.bcf &
bcftools view -S filtering_lists/lib_three_ids.filtered.txt -Ob -o lib_three.invariant.bcf filtering_intermediates_invariant/SNPs_INFO_and_missing_filtered.invariant_sites.bcf &
bcftools view -S filtering_lists/lib_four_ids.filtered.txt -Ob -o lib_four.invariant.bcf filtering_intermediates_invariant/SNPs_INFO_and_missing_filtered.invariant_sites.bcf &

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
bcftools view all_libs.invariant.DP-GQ_filtered.bcf | grep -m 1 "^#CHROM" | grep -o "NG" | wc -l

#filter sites on missingness > 25%
bcftools view -e 'F_MISSING>0.25' -Ob -o FINAL_invariant.bcf all_libs.invariant.DP-GQ_filtered.bcf

#check number of sites removed (need to run when on power)
bcftools view FINAL_invariant.bcf | grep -v "^#" | wc -l
#41026896
bcftools view all_libs.invariant.DP-GQ_filtered.bcf | grep -v "^#" | wc -l
# should be 42335985 methinks; it is
# 42335985 - 41026896 = 1309089
#


#count hits to mt genome
bcftools view FINAL_invariant.bcf | grep "^tig00000405_pilon"  | wc -l
#8338 remove these
bcftools view FINAL_invariant.bcf | grep "^tig00000405_pilon"  | grep ":RGQ" | wc -l
# 8262
bcftools view FINAL_invariant.bcf | grep "^tig00000405_pilon"  | grep ":GQ" | wc -l
# 76

bcftools index FINAL_invariant.bcf
bcftools view FINAL_invariant.bcf | grep -v "^tig00000405_pilon" | bcftools view -Ob -o FINAL_invariant.nuclear.bcf

bcftools view FINAL_invariant.nuclear.bcf | grep "^tig00000405_pilon"  | wc -l
bcftools view FINAL_invariant.nuclear.bcf | grep -m 1 "^#CHROM" | grep -o "NG" | wc -l
#117
bcftools view FINAL_invariant.bcf | grep -m 1 "^#CHROM" | grep -o "NG" | wc -l
#117

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
grep "^tig00000405_pilon" FINAL_snp.vcf | wc -l
#76 remove these
grep -v "^tig00000405_pilon" FINAL_snp.vcf > FINAL_snp.nuclear.vcf
bcftools view -Ob -o FINAL_snp.nuclear.bcf FINAL_snp.nuclear.vcf

#there are some samples of different individuals from the same tree;
# N149, 118 (ANF1 1); NG114, 152 (ANF1 10);
#look at which ones have better completeness
cd ~/repo/neonectria_SNP/data/Nf/final_tables/
grep -m 1 "^#CHROM" FINAL_snp.nuclear.vcf | grep -o "NG" | wc -l
vcftools --vcf FINAL_snp.nuclear.vcf --missing-indv
#117

grep -E 'NG149|NG118|NG114|NG152' out.imiss

#remove dups from vcfs as well
printf 'NG149\nNG152' > ../filtering_lists/dups.txt
mkdir rm_dups

less ../filtering_lists/dups.txt

bcftools view -S ^../filtering_lists/dups.txt -Oz -o rm_dups/FINAL_invariant.nuclear.vcf.gz FINAL_invariant.nuclear.bcf &
bcftools view -S ^../filtering_lists/dups.txt -Oz -o rm_dups/FINAL_snp.IBD_analyses.vcf.gz FINAL_snp.nuclear.bcf &

gunzip -c rm_dups/FINAL_invariant.nuclear.vcf.gz | grep -m 1  "^#CHROM" | grep -o "NG" | wc -l
#115
gunzip -c rm_dups/FINAL_snp.IBD_analyses.vcf.gz | grep -m 1  "^#CHROM" | grep -o "NG" | wc -l
#115
bcftools view rm_dups/FINAL_invariant.nuclear.vcf.gz | grep -v '^#' | wc -l
# 41040857
bcftools view rm_dups/FINAL_snp.IBD_analyses.vcf.gz | grep -v '^#' | wc -l
#1116141
bcftools view rm_dups/FINAL_snp.IBD_analyses.vcf.gz > rm_dups/FINAL_snp.IBD_analyses.vcf
python ~/repo/vcf2phylip/vcf2phylip.py --input rm_dups/FINAL_snp.IBD_analyses.vcf -m 1
#Total of genotypes processed: 1116141
#Genotypes excluded because they exceeded the amount of missing data allowed: 0
#Genotypes that passed missing data filter but were excluded for being MNPs (multiple nuc poly): 119050
#SNPs that passed the filters: 997091
rm rm_dups/FINAL_snp.IBD_analyses.vcf

#rm singletons
cd rm_dups
vcftools --gzvcf FINAL_snp.IBD_analyses.vcf.gz --mac 2 --recode --out FINAL_snp.mac_ge2
#kept 447150 out of a possible 1116141 Sites
bgzip FINAL_snp.mac_ge2.recode.vcf

#bialleles
bcftools view -m2 -M2 -v snps FINAL_snp.IBD_analyses.vcf.gz -Oz -o FINAL_snp.biallele.vcf.gz
#convert to diploid for calculation of tajD with vcf-kit
# note the -- indicates start of plugin options
bcftools +fixploidy -Oz -o FINAL_snp.biallele.TajDdiploid.vcf.gz FINAL_snp.biallele.vcf.gz -- -f 2
########
# the snpR package is absolute garbage, don't bother

#remove first space in header
bgzip FINAL_snp.biallele.snpR.vcf
gunzip FINAL_snp.biallele.snpR.vcf

gunzip -c FINAL_snp.biallele.vcf.gz | grep -v "^#" | cut -f-5,7- > FINAL_snp.biallele.no_qual.data.only
gunzip -c FINAL_snp.biallele.vcf.gz | grep "^#" > FINAL_snp.biallele.no_qual.header
#manually remove the QUAL from the header
cat FINAL_snp.biallele.no_qual.header FINAL_snp.biallele.no_qual.data.only > FINAL_snp.biallele.no_qual.vcf
bgzip FINAL_snp.biallele.no_qual.vcf

#
#LD filter
bcftools +prune -m 0.5 -w 10000 -Oz -o FINAL_snp.mac_ge2.LD.pca_analyses.vcf.gz FINAL_snp.mac_ge2.recode.vcf.gz &

#
#bialleles & LD
bcftools view -m2 -M2 -v snps FINAL_snp.mac_ge2.recode.vcf.gz | bcftools +prune -m 0.5 -w 10000 -Oz -o FINAL_snp.mac_ge2.biallele.LD.structure_analyses.vcf.gz &

#
#bialleles (on MAC filtered)
bcftools view -m2 -M2 -v snps FINAL_snp.mac_ge2.recode.vcf.gz -Oz -o FINAL_snp.mac_ge2.biallele.gwas_analyses.vcf.gz &

########################
# After final filtration we filter the metadata to only the samples retained
#
# metadata_collate/filter_to_retained_samples.R

#############################################################################
#############################################################################
# After performing pop gen analysis
# we discover that most of the samples are in a single PCA/STRUCTURE/ADMIXTURE cluster
# we will isolate these samples for LFMM analyses (and pixy genome scan with invariant)
# to avoid issues with "allele surfing"

#

cd ~/repo/neonectria_SNP/data/Nf/final_tables/rm_dups/
bcftools view -S ../../filtering_lists/PCA_cluster_2_ids.txt FINAL_snp.mac_ge2.biallele.gwas_analyses.vcf.gz -Oz -o FINAL_snp.gwas_analyses.cluster2.vcf.gz
bcftools view -S ../../filtering_lists/PCA_cluster_2_ids.txt FINAL_invariant.nuclear.vcf.gz -Oz -o FINAL_invariant.nuclear.cluster2.vcf.gz

cd ~/repo/neonectria_SNP/data/Nf/final_tables/rm_dups/
bcftools view -S ../../filtering_lists/PCA_cluster_2_ids_no_NC.txt FINAL_snp.mac_ge2.biallele.gwas_analyses.vcf.gz -Oz -o FINAL_snp.gwas_analyses.cluster2_no_NC.vcf.gz
bcftools view -S ../../filtering_lists/PCA_cluster_2_ids_no_NC.txt FINAL_invariant.nuclear.vcf.gz -Oz -o FINAL_invariant.nuclear.cluster2_no_NC.vcf.gz
