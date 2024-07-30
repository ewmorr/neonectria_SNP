require(dplyr)
require(ggplot2)
require(stringr)

source("~/repo/neonectria_genome_reseq_10072020/R_scripts/ggplot_theme.txt")

r_tab = read.table("~/Nf_SPANDx_all_seqs/Outputs/Master_vcf/out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.LD_decay.ld.gz", header = T)

#remove uneeded cols to save mem
r_tab = subset(r_tab, select = -c(SNP_A, CHR_B, SNP_B))
#calculate SNP dist
r_tab$distance = r_tab$BP_B - r_tab$BP_A

###############
#calculate means
#different dists at 1K, 2.7.5K, 7.5K, 10K
###############

#############
#1KB mean
#############
print("starting 1K means")
summarize_dist = 1000
r_tab$distc <- cut(
  r_tab$distance,
  breaks=seq(
    from=0,
    to=ceiling(max(r_tab$distance)/summarize_dist)*summarize_dist,
    by=summarize_dist
    ) 
  )

r_tab_avg <- r_tab %>% group_by(distc, CHR_A) %>% summarise(mean=mean(R2),median=median(R2))

#calc mean distance
r_tab_avg <- r_tab_avg %>% mutate(
    start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
    end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
    mid=start+((end-start)/2)
    )

write.table(r_tab_avg, "LD.1KB_avg.txt", sep = "\t", quote = F, col.names = T, row.names = F)

