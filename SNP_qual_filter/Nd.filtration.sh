
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


######################################

######################################
######################################
######################################
#Invariant sites filtration
#
# 1. GATK filters
# 2. DP (based on variant site means per lib), GQ, abd RGQ filters
# 3. remove lists of indv and sites
# 
# The invariants site vcf is on the server...

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
