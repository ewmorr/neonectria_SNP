
#run locally
conda activate bcftools

cd ~/repo/neonectria_SNP/data/shared_buscos/

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
# kept 147 of 147 indv
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

#vcftools --vcf INFOfilters.removed.recode.vcf --out lib_one --recode --keep lib_one_ids.txt
#vcftools --vcf INFOfilters.removed.recode.vcf --out lib_two --recode --keep lib_two_ids.txt
#vcftools --vcf INFOfilters.removed.recode.vcf --out lib_three --recode --keep lib_three_ids.txt
#vcftools --vcf INFOfilters.removed.recode.vcf --out lib_four --recode --keep lib_four_ids.txt
#vcftools --vcf INFOfilters.removed.recode.vcf --out lib_five --recode --keep lib_five_ids.txt

#Calculate DP values in R before further processing

#- locus DP across libraries (these have been updated to shared_buscos from the Nf copy over)
## from (Nd_variant_quality_exploratory.r)
#
#    libOne  libTwo  libThree    libFour libFive
#max (min is 0)   99.71014 164.9512  191.4286   92.68   83.2
#mean    8.1176    23.70867   45.57019  38.26835    23.46229
#median  7.869565    23.90244   47.71429  41.24 23.6
#sd  1.993602 4.338694    10.62207    9.955711  4.583732
#mean + 2*sd    12.1048    32.38606    66.81434   58.17977  32.62976
#mean - 2*sd 1.852007   15.03128 24.32604   18.35692    14.29483

########################################################################
########################################################################
#NOTE: the point of the max DP filter is to pick up repeat regions. These BUSCOs are all 
# single copy orthologs so it's not needed. can actually fully skip the lib-based
# depth filter and just apply the min DP filter to the full lib. (that's also
# presumably why the SD is so tight and we filter out a lot of sites when we use mean+2*SD)
# 
#let's go with a minimum DP of 2 and minimum GQ of 30. Because we are haploid we have greater confidence as well... i.e., we do not need to predict/worry about heterozygotes (Lou et al. 2021)
#We will also pick up low confidence SNPs with percent NA filters (Lou et al. 2021)


#vcftools --vcf lib_one.recode.vcf --out lib_one.DPGQ_filter --recode --minDP 2 --minGQ 30 --maxDP 12
#vcftools --vcf lib_two.recode.vcf --out lib_two.DPGQ_filter --recode --minDP 2 --minGQ 30 --maxDP 32
#vcftools --vcf lib_three.recode.vcf --out lib_three.DPGQ_filter --recode --minDP 2 --minGQ 30 --maxDP 66
#vcftools --vcf lib_four.recode.vcf --out lib_four.DPGQ_filter --recode --minDP 2 --minGQ 30 --maxDP 58
#vcftools --vcf lib_five.recode.vcf --out lib_five.DPGQ_filter --recode --minDP 2 --minGQ 30 --maxDP 32

vcftools --vcf INFOfilters.removed.recode.vcf --out INFOfilters.removed.DPGQ_filter --recode --minDP 2 --minGQ 30 


############################
#below needs to be updated for this dataset
#count number flagged missing
grep -v "^#" INFOfilters.removed.DPGQ_filter.recode.vcf | wc -l
#755215 * 147 = 111016605 (total sites)
grep -v "^#" INFOfilters.removed.recode.vcf | grep -o "\s\.:" | wc -l
#3742424 (prefilter)
grep -v "^#" INFOfilters.removed.DPGQ_filter.recode.vcf | grep -o "\s\.\/\.:" | wc -l
#4063035 (post filter)
#4063035-3742424 = 320611
# 320611/111016605 = 0.002887955
# 4063035/111016605 = 0.03659844



#missing data filters after marking sites missing based on depth
#We already excluded all of the individuals from SNP calling
#that were excluded from the full datasets
# so we will simply filter at 25% missing for SNPs
# Still check ind missingness here

# the LYLF paper applies: missing inds per genotype (geno) < 50%, missing data 
# individuals (imiss) < 90%; geno < 40%, imiss < 70%; geno < 30%, imiss < 50%; 
# geno < 5%, imiss < 25%


#get percent missing per site and indv to look at distribution to identify
# levels for filtering
# subbing . back in for ./. before calculate
sed 's@\./\.@\.@g' INFOfilters.removed.DPGQ_filter.recode.vcf > INFOfilters.removed.DPGQ_filter.hap.vcf

vcftools --vcf INFOfilters.removed.DPGQ_filter.hap.vcf --missing-indv & 
vcftools --vcf INFOfilters.removed.DPGQ_filter.hap.vcf --missing-site & 

mv out.imiss INFOfilters.removed.DPGQ_filter.imiss 
mv out.lmiss INFOfilters.removed.DPGQ_filter.lmiss 

###########################################
#First pass -- all individuals are below 15% 
# apply a single hard filter to the SNPs to 25% 

#somehow a comma jumped back in here. It's in the SF INFO field and looks like
#must have been inserted by vcftools (eyeroll)
perl -pe 's/(".+?[^\\]")/($ret = $1) =~ (s#,##g); $ret/ge' INFOfilters.removed.DPGQ_filter.hap.vcf > INFOfilters.removed.DPGQ_filter.hap.noComma.vcf

lmiss=0.25

awk -v var="$lmiss" '$6 > var' INFOfilters.removed.DPGQ_filter.lmiss | cut -f1,2 > NA.${lmiss}.sites
tail -n +2 NA.${lmiss}.sites > temp.sites && mv temp.sites NA.${lmiss}.sites
wc -l NA.${lmiss}.sites
wc -l INFOfilters.removed.DPGQ_filter.lmiss
#755216-12657 = 742559
#12657/755216 = 0.01675944
#removing loci with >25% missing sites (70947 loci; 488725 total; 14.5%)
vcftools --vcf INFOfilters.removed.DPGQ_filter.hap.noComma.vcf --exclude-positions NA.${lmiss}.sites --recode --out NA_filter.loc_gt_$lmiss
# kept 742558 out of a possible 755215 Sites

#Recalculate missingness
vcftools --vcf NA_filter.loc_gt_$lmiss.recode.vcf --missing-indv &
vcftools --vcf NA_filter.loc_gt_$lmiss.recode.vcf --missing-site &

mv out.imiss NA_filter.loc_gt_$lmiss.imiss
mv out.lmiss NA_filter.loc_gt_$lmiss.lmiss


##########################
#cleanup and archive
#
#Final table (pass to filtering for polyalleles, MAC, LD as necessary depending on analysis)
cd ~/repo/neonectria_SNP/data/shared_buscos
cp NA_filter.loc_gt_0.25.recode.vcf FINAL_snp.vcf
bcftools view -Ob -o FINAL_snp.bcf FINAL_snp.vcf &

#move intermediate files and zip
mkdir filtering_intermediates
mv NA* filtering_intermediates
mv INFOfilters* filtering_intermediates

tar -czvf filtering_intermediates.tar.gz filtering_intermediates &

######################################
######################################

#Once the filtering is complete we can remove the filtering_intermediates folder
rm -R filtering_intermediates

mkdir unfiltered_vcfs
mv out* unfiltered_vcfs
bgzip unfiltered_vcfs/out.vcf
tar -cvf unfiltered_vcfs.tar unfiltered_vcfs
rm -R unfiltered_vcfs

#clean up final tables
mkdir final_tables 
mv FINAL_snp* final_tables


#######################
#######################

#final filters for analyses
cd final_tables

#there are some Nd samples of different individuals from the same tree or the same indv sequenced twice;
# NG163, NG20 (MI1 1.1.1); NG27, NG144 (MI1 2);
#look at which ones have better completeness
cd ~/repo/neonectria_SNP/data/shared_buscos/final_tables/
vcftools --vcf FINAL_snp.vcf --missing-indv
grep -E 'NG163|NG20|NG27|NG144' out.imiss

#remove dups from vcfs as well (these are same samples as full dataset)

#there are some samples of different individuals from the same tree;
# N149, 118 (ANF1 1); NG114, 152 (ANF1 10);
#look at which ones have better completeness
grep -E 'NG149|NG118|NG114|NG152' out.imiss
#the missingness is low in all cases so we rm the same samples as the full dataset (would be 118 instead of 149)
printf 'NG20\nNG27\nNG149\nNG152' > ../filtering_lists/dups.txt

mkdir rm_dups

less ../filtering_lists/dups.txt

bcftools view -S ^../filtering_lists/dups.txt -Oz -o rm_dups/FINAL_snp.IBD_analyses.vcf.gz FINAL_snp.bcf 

bcftools view rm_dups/FINAL_snp.IBD_analyses.vcf.gz | grep -v '^#' | wc -l
#742558

########################
# After final filtration the Nf/Nd metadata are filtered to only the samples retained
# metadata_collate/filter_to_retained_samples.R
# 
# For the full set can combine the species metadata as the same samples are present
