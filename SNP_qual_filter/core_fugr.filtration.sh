
#run locally
conda activate bcftools

cd ~/repo/neonectria_SNP/data/Fugr1_ref/

QDparam=2
MQparam=40
FSparam=60

gatk CreateSequenceDictionary -R core.Fusgr1_AssemblyScaffolds_Repeatmasked-N.fasta
samtools faidx core.Fusgr1_AssemblyScaffolds_Repeatmasked-N.fasta

gatk VariantFiltration -R core.Fusgr1_AssemblyScaffolds_Repeatmasked-N.fasta -V out.vcf -O INFOfilters.vcf \
    --filter-name "MQ40" -filter "MQ < 40.0" \
    --filter-name "QD2" -filter "QD < 2.0" \
    --filter-name "FS60" -filter "FS > 60.0"
    

#remove commas from vcf INFO field descriptions before proceeding with vcftools
# https://unix.stackexchange.com/questions/48672/only-remove-commas-embedded-within-quotes-in-a-comma-delimited-file

perl -pe 's/(".+?[^\\]")/($ret = $1) =~ (s#,##g); $ret/ge' INFOfilters.vcf > INFOfilters.noComma.vcf


# For the diversity analyses we will keep all mac and all allele counts. We have other qual fiulters to depend on to filter out bad calls, and these types of calls are obviously important to calculation of pi and Tajima's D, for example. Also see Lou et al 2021 which points out these should not be removed for div analyses

vcftools --vcf INFOfilters.noComma.vcf --out INFOfilters.removed --recode --remove-filtered-all
# kept 147 of 147 indv
# After filtering, kept 52383 out of a possible 53288 Sites
# 53288 - 53283 = 5

#there are some samples of different individuals from the same tree;
# Nf N149, 118 (ANF1 1); NG114, 152 (ANF1 10);
# Nd NG163, NG20 (MI1 1.1.1); NG27, NG144 (MI1 2);
#rm the same samples as the full dataset 
mkdir filtering_lists
printf 'NG20\nNG27\nNG149\nNG152' > filtering_lists/dups.txt
bcftools view -S ^filtering_lists/dups.txt -Ov -o INFOfilters.removed.rm_dups.vcf INFOfilters.removed.recode.vcf 

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

vcftools --vcf INFOfilters.removed.recode.vcf --out lib_one --recode --keep lib_one_ids.txt &
vcftools --vcf INFOfilters.removed.recode.vcf --out lib_two --recode --keep lib_two_ids.txt &
vcftools --vcf INFOfilters.removed.recode.vcf --out lib_three --recode --keep lib_three_ids.txt &
vcftools --vcf INFOfilters.removed.recode.vcf --out lib_four --recode --keep lib_four_ids.txt &
vcftools --vcf INFOfilters.removed.recode.vcf --out lib_five --recode --keep lib_five_ids.txt &

#Calculate DP values in R before further processing

#- locus DP across libraries (these have been updated to shared_buscos from the Nf copy over)
## from (Nd_variant_quality_exploratory.r)
#
#    libOne  libTwo  libThree    libFour libFive
#max (min is 0)   647.8116 1319.537  2678.571   1243.96   865.6
#mean    8.658848    21.37757   37.42161  29.77809    19.59264
#median  9.565217    23.41463   38.42857  30.96 21
#sd  7.185259 15.74534    28.72833    18.44462  10.81086
#mean + 2*sd    23.02937    52.86825    94.87828   66.66734  41.21436
#mean - 2*sd all negative

########################################################################
# 
#let's go with a minimum DP of 2 and minimum GQ of 30. Because we are haploid we have greater confidence as well... i.e., we do not need to predict/worry about heterozygotes (Lou et al. 2021)
#We will also pick up low confidence SNPs with percent NA filters (Lou et al. 2021)


vcftools --vcf lib_one.recode.vcf --out lib_one.DPGQ_filter --recode --minDP 2 --minGQ 30 --maxDP 23 &
vcftools --vcf lib_two.recode.vcf --out lib_two.DPGQ_filter --recode --minDP 2 --minGQ 30 --maxDP 52 &
vcftools --vcf lib_three.recode.vcf --out lib_three.DPGQ_filter --recode --minDP 2 --minGQ 30 --maxDP 94 &
vcftools --vcf lib_four.recode.vcf --out lib_four.DPGQ_filter --recode --minDP 2 --minGQ 30 --maxDP 66 &
vcftools --vcf lib_five.recode.vcf --out lib_five.DPGQ_filter --recode --minDP 2 --minGQ 30 --maxDP 41 &

#vcftools --vcf INFOfilters.removed.recode.vcf --out INFOfilters.removed.DPGQ_filter --recode --minDP 2 --minGQ 30 


############################
#below needs to be updated for this dataset
#count number flagged missing
grep -v "^#" lib_one.recode.vcf | wc -l
wc -l lib_one_ids.txt
#52383 * 99 = 5185917 (total sites)
grep -v "^#" lib_one.recode.vcf | grep -o "\s\.:" | wc -l
#741498 (prefilter)
zgrep -v "^#" lib_one.DPGQ_filter.recode.vcf.gz | grep -o "\s\.\/\.:" | wc -l
#895458 (post filter)v
#895458-741498 = 153960
# 153960/5185917 = 0.0296881
# 895458/5185917 = 0.1726711

# merge the vcfs from the four libs.
# need to bgzip and tabix first

for i in one two three four five
do(
    bgzip lib_${i}.DPGQ_filter.recode.vcf
    tabix lib_${i}.DPGQ_filter.recode.vcf.gz
)
done

vcf-merge lib_one.DPGQ_filter.recode.vcf.gz \
    lib_two.DPGQ_filter.recode.vcf.gz \
    lib_three.DPGQ_filter.recode.vcf.gz \
    lib_four.DPGQ_filter.recode.vcf.gz \
    lib_five.DPGQ_filter.recode.vcf.gz > all_libs.DPGQ_filter.recode.vcf


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
sed 's@\./\.@\.@g' all_libs.DPGQ_filter.recode.vcf > INFOfilters.removed.DPGQ_filter.hap.vcf

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

# the data is high missingness. Start with 50% filter. This may be all we do
# # It's likely due to mapping rates against FUGR?
lmiss=0.50

awk -v var="$lmiss" '$6 > var' INFOfilters.removed.DPGQ_filter.lmiss | cut -f1,2 > NA.${lmiss}.sites
tail -n +2 NA.${lmiss}.sites > temp.sites && mv temp.sites NA.${lmiss}.sites
wc -l NA.${lmiss}.sites
wc -l INFOfilters.removed.DPGQ_filter.lmiss
#52384-9843 = 42541
#9843/52384 = 0.1879009
#removing loci with >50% missing sites (70947 loci; 488725 total; 14.5%)
vcftools --vcf INFOfilters.removed.DPGQ_filter.hap.noComma.vcf --exclude-positions NA.${lmiss}.sites --recode --out NA_filter.loc_gt_$lmiss
# kept 42540 out of a possible 52383 Sites

#Recalculate missingness
vcftools --vcf NA_filter.loc_gt_$lmiss.recode.vcf --missing-indv &
vcftools --vcf NA_filter.loc_gt_$lmiss.recode.vcf --missing-site &

mv out.imiss NA_filter.loc_gt_$lmiss.imiss
mv out.lmiss NA_filter.loc_gt_$lmiss.lmiss


##########################
#cleanup and archive
#
#Final table (pass to filtering for polyalleles, MAC, LD as necessary depending on analysis)
cd ~/repo/neonectria_SNP/data/Fugr1_ref
cp NA_filter.loc_gt_0.50.recode.vcf FINAL_snp.vcf
bcftools view -Ob -o FINAL_snp.bcf FINAL_snp.vcf &

#move intermediate files and zip
mkdir filtering_intermediates
mv NA* filtering_intermediates
mv INFOfilters* filtering_intermediates
mv lib* filtering_intermediates
mv all_libs* filtering_intermediates
rm filtering_intermediates/*log

tar -czvf filtering_intermediates.tar.gz filtering_intermediates 

#get the INFO filtered sites list before moving on
grep -v "^#" filtering_intermediates/INFOfilters.noComma.vcf | grep -v "PASS" | cut -f 1,2 > filtering_lists/INFO_filters.sites
grep -v "^#" filtering_intermediates/INFOfilters.noComma.vcf | grep "PASS" | cut -f 1,2 | wc -l
grep -v "^#" filtering_intermediates/INFOfilters.removed.recode.vcf | grep "PASS" | cut -f 1,2 | wc -l

cp filtering_intermediates/NA.0.50.sites filtering_lists
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
# Invariant site filters
#

cat filtering_lists/NA.0.50.sites filtering_lists/INFO_filters.sites > filtering_lists/all_to_filter.sites
tar -zxvf unfiltered_vcfs.tar unfiltered_vcfs/out.invariant_sites.vcf.gz
bcftools view -S ^filtering_lists/dups.txt -T ^filtering_lists/all_to_filter.sites -Oz -o SNPs_INFO_and_missing_filtered.invariant_sites.vcf.gz unfiltered_vcfs/out.invariant_sites.vcf.gz 

gunzip -c SNPs_INFO_and_missing_filtered.invariant_sites.vcf.gz | grep -v "^#" | wc -l
#324121
gunzip -c unfiltered_vcfs/out.invariant_sites.vcf.gz | grep -v "^#" | wc -l
# 334869 - 324121 = 157997

bcftools view SNPs_INFO_and_missing_filtered.invariant_sites.vcf.gz -Ob -o SNPs_INFO_and_missing_filtered.invariant_sites.bcf 

#filtering for invariant sites
bcftools +setGT SNPs_INFO_and_missing_filtered.invariant_sites.bcf  -Ob -o invariant.DP-GQ_filter.bcf -- -t q -n . -i 'FMT/DP<2 | FMT/DP>46 | GQ<30 | RGQ<30' &
# Filled 6463910 alleles

#filter sites on missingness > 50%
bcftools view -e 'F_MISSING>0.5' -Ob -o FINAL_invariant.bcf invariant.DP-GQ_filter.bcf

#check number of sites removed 
bcftools view -Ov FINAL_invariant.bcf | grep -v "^#" | wc -l
#302781
bcftools view -Ov invariant.DP-GQ_filter.bcf | grep -v "^#" | wc -l
# should be 324121 methinks; it is
tar -xvf filtering_intermediates.tar.gz
mv SNPS* filtering_intermediates
mv invariant* filtering_intermediates
tar -cvf filtering_intermediates.tar.gz filtering_intermediates
rm -R filtering_intermediates

