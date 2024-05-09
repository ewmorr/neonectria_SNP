# Repo for analysis of genome resequencing data of Neonectria faginata and Neonectria ditissima
#### Workflows performed on UNH Premise or locally using R as appropriate

## Initial processing and SNP calling performed as described in `README.SPANDx_SNP_calling.md`

## SNP and sample filtering

### Filter out SNPs that SPANDx FAIL filtering
```
cd ~/Nf_SPANDx_all_seqs/Outputs/Master_vcf
ls ../.. | grep fastq | wc -l #250/2 = 125
grep "##\|#\|PASS" out.filtered.vcf > out.filtered.PASS.vcf
grep -v "##\|#" out.filtered.PASS.vcf | wc -l
```
530256 SNPs remaining x 125 samples = 66282000
```
grep -o "\s\.:" out.filtered.PASS.vcf | wc -l
```
2920386 NA sites
2920386/66282000 = 4.41%

#### Note that SPANDx performs variant level filtering (i.e., across samples) based on depth calibrated quality score (QD, QDFilter), root mean square mapping quality (MQ) and Fisher Strand bias (Fisher's exact test on reads supporting different base call on fwd vs rev strand [high FS is worse], QD > 10 (GATK rec is 2); MQ > 30 (GATK rec is 40); FS < 60 (GATK rec is 60). Note that QD only applies to variant calls (so does not affect invariant sites)
#### See here for variant filtering reccomendations from Broad institute https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
#### We need to examine the effect of the SPANDx filtering and possibly perform a different filtering based on the GATK recs
- the largest effect will be on the change in QD,which retains more SNPs with the GATK rec. There is a long tail (skew) and the SPANDx filter of QD > 10 is about halfway through. Not a huge effect (visually) and we can likely take the GATK rec, based on increasing the MQ filter (note that a QUAL filter of min 30 has been performed automatically as GATK rec'd; this is the locus level quality score, so feeds QD), and then performing subsequent depth based filtering. Filter at GATK recs across the board

Terminology from "loci you're looking for" paper
locus, SNP (a position in the genome across individuals)
individual (a sample)
genotype (ind x locus)

1. First filter loci based on GATK recs, minor allele count (mac >= 2), biallelic only (the latter two with vcftools) AND
2. apply depth filters either across sequencing libraries or per individual; can apply DP filters based on individual mean or library wide mean; may make sense to filter whole loci based on the lib-wise mean, and filter individual genotypes based on the individual mean
- There are some strong library based effects on both locus DP and individual DP. Follow above, i.e., filter whole loci based on lib-wise 
- DP across libraries
```
    libOne  libTwo  libThree
range   0.0000-840.7727 0.000-1902.646  0.000-2264.636
mean    8.09    24.48   29.08
median  6.54    23.77   28.81

```
--DP genotype DP filters should be applied on VCF that is split by library (as there are library depth effects), and then the split VCF recombined before performing the following missingness filtering (vcf-merge to cat samples with the same loci, vcf-concat to cat loci with the same samples [e.g., different chromosomes], must be bgzipped and tabixed)
3. apply iterative missing data filters, e.g., the LYLF paper applies: geno < 50%, imiss < 90%; geno < 40%, imiss < 70%; geno < 30%, imiss < 50%; geno < 5%, imiss < 25%
- note LYLF paper applies depth and qual based filters at genotype level (minDP 5, minQUAL 20), applies mean DP (15) and mac (3), THEN filters by iterative NA filters, then appliers INFO based filters (i.e., the GATK recs), then applies the final missingness filter. WHY? It seems that GATK hard-filtering should be applied first as these are the recommended filters for low-qual SNPs and can therefore affect missingness at both locus and individual


## For invariant calls set
### Invariant calls do not get QUAL scores (reserved for variant sites) but instead get variant level (i.e., genotype x ind) GT:DP:RGQ. RGQ is the reference genotype quality, or the Phred scaled probability the call is correct. We can therefore filter individual variants based on GQ. 
### RGQ = -10 x log10(p). So that RGQ = 30 is a 1 percent prob the call is *incorrect*, p = 0.001 == RGQ = 30, p = 0.0001 == RGQ = 40,
### Note that this is synonymous to the site level QUAL (but at the sample level) so scaling GQ by per sample DP makes some sense, but the threshold of 2 does not seem to work well...

#### For testing of RGQ and ref genotype calls we want to compare DP RGQ and whether a GT is called (i.e., not NA)
### GT:DP:RGQ
#### Can plot this against calls
```
    #rm notes/comments
    #get only non variants 
    #the last grep needs -E (extended regexp) to correctly interpret the `+` wildcard with \d
    
gunzip -c out.invariant_sites.vcf.gz |
    grep -v ^# | 
    grep RGQ | 
    grep -o -E '\s\.:\d+:\d+' | 
    wc -l
#148394771

gunzip -c out.invariant_sites.vcf.gz |
    grep -v ^# | 
    grep RGQ | 
    grep -o -E '\s\.:\d+:\d+' | 
    cut -d ':' -f 2 | sort | uniq
#NA site depth ranges from 0 to 2915, worthwhile to plot
#maybe plot DP versus RGQ and color by NA or "called" (or maybe REF and ALT)

#write out files for plotting. probably most efficient to split the metrics outside of R and then pull in and join to plot
gunzip -c out.invariant_sites.vcf.gz |
    grep -v ^# | 
    grep RGQ | 
    grep -o -E -e '\s(\.|0):\d+:\d+' | #note there will be no ALT calls bc these are invariant sites (duh)
    cut -d ':' -f 1 > NAvREF.invariant_sites.genoXindv.txt
    
    
gunzip -c out.invariant_sites.vcf.gz |
    grep -v ^# | 
    grep RGQ | 
    grep -o -E -e '\s(\.|0):\d+:\d+' | #note there will be no ALT calls bc these are invariant sites (duh)
    cut -d ':' -f 2 > DP.invariant_sites.genoXindv.txt
    
gunzip -c out.invariant_sites.vcf.gz |
    grep -v ^# | 
    grep RGQ | 
    grep -o -E -e '\s(\.|0):\d+:\d+' | #note there will be no ALT calls bc these are invariant sites (duh)
    cut -d ':' -f 3 > RGQ.invariant_sites.genoXindv.txt
```
These files are big (like 18G). Will need to upload to premise for plotting. Also 5152446314 lines. Probably need to subsample. use gshuf (mac) or shuf (linux)) OR Maybe just start with one scf. 560M lines would be about 100 samples of the first scf
```
gunzip -c out.invariant_sites.vcf.gz |
    grep -v ^# | 
    grep RGQ | 
    grep -o -E -e '\s(\.|0):\d+:\d+' |
    shuf -n 10000000 > NA_DP_RGQ.invariant_sites.genoXindv.subsample.txt

cut -d ':' -f 1 NA_DP_RGQ.invariant_sites.genoXindv.subsample.txt > NAvREF.invariant_sites.genoXindv.subsample.txt
cut -d ':' -f 2 NA_DP_RGQ.invariant_sites.genoXindv.subsample.txt > DP.invariant_sites.genoXindv.subsample.txt
cut -d ':' -f 3 NA_DP_RGQ.invariant_sites.genoXindv.subsample.txt > RGQ.invariant_sites.genoXindv.subsample.txt
```


#### We also want to plot depth across the entire call set
#### we need to replace dot's in the info call with DP=0 to process correctly
#### first cut out the INFO field. The ref only calls will only have DP (we think) or dot (which should be replaced with zero with sed), then can split out the rest of the DP calls with SED
#### Also want to pull out the position of the base call for plotting (and probably the chromosome)
```
gunzip -c out.invariant_sites.vcf.gz |
    grep -v ^# | 
    cut -f 8 | #is this the right field
    sed  's/^\.$/DP=0/' | #replace dot with DP=0
    grep  -o -E 'DP=\d+' | #grep pull out the DP from multi field INFO lines
    grep -o -E '\d+' > DP.invariant_sites.txt #extract the numeric only
    
gunzip -c out.invariant_sites.vcf.gz |
    grep -v ^# | 
    cut -f 1-2 > scf-pos.invariant_sites.txt
```

####################
# Combined call set filtering
The INFO fields for both variant and invariant can be filtered as normal because the scoring categories are unique to variant/invaraiant. However, due to the iterative missing data filters the locus and individual calls will both (presumably) be affected by including both variant and invariant sites. We perform variant filtering as above, and then note the individuals and loci that are exlcuded from the final table. We then apply name based filters as opposed to filtering based on calculated missingness.


```
#iterate on the cutoff vars (lmiss == locus; imiss == individual)
lmiss=
imiss=

vcftools -vcf target.vcf --missing-site
awk '$5 > $lmiss' out.lmiss | cut -f1 > NA.${lmiss}.sites
vcftools -vcf target.vcf --exclude NA.${lmiss}.sites --recode --out target.loc${lmiss}.vcf

vcftools -vcf target.site${missingnessCutoff}.vcf --missing-indv
awk '$5 > $imiss' out.imiss | cut -f1 > NA.${imiss}.indv
vcftools -vcf target.vcf --remove NA.${imiss}.sites --recode --out target.loc-$lmiss.ind-${imiss}.vcf

```
The above run with repeated calls in the iterative filtering along with individual NA filtering, and then cat the `NA.*.sites` and `NA.*.indv` files and run and --exclude and a --remove on the combined variant/invariant file

*Note that vcftools removes data, NOT flags data as not PASS*




