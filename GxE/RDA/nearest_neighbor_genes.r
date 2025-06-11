library(dplyr)
library(data.table)

#############################
#FUNS
#############################

#function for nearest neighbor
#we will have already searched for SNPs *within* genes so the we can then search for start sites > then pos else stop sites < pos for the nearest neighbor match
#genes down stream of SNPs will have positive distance and upstream will have negative distance
#i.e., start - SNPpos OR stop - SNPpos
nearest_gene = function(snpPos = NULL, gff = NULL){
    #new col for distances
    gff$dist = vector(mode = "numeric", length = nrow(gff))
    #calculate distances
    gff$dist = ifelse(gff$start > snpPos, gff$start - snpPos, gff$stop - snpPos)
    #filter to minimum absolute distance and return the resulting df
    #dplyr way
    #return( filter(gff, abs(dist) == min(abs(gff$dist))) )
    #DT way
    return(
        gff[
            abs(dist) == min(abs(gff$dist))
        ]
    )
}

############################################


##########################
#Import data
##########################

#import gff (gene positions) contig name, gene ID, start site, stop site
gff = read.csv("data/Nf/funann.mRNA.gff3", sep = "\t")
head(gff)

#import SNP pos results
#
set_me = "no_NC.reduced"

sig_SNPs = read.csv("data/Nf/GxE/RDA/no_NC.rdadapt_reduced_env.csv")
colnames(sig_SNPs)
nrow(sig_SNPs)
sig_SNPs$scaffold = sub("_pilon", "", sig_SNPs$scaffold)
head(sig_SNPs)
#filter for SNPs with sig relationship to variable
sig_SNPs.sig = sig_SNPs %>% filter(rdadapt_sig == "sig") %>% select(scaffold, position)
nrow(sig_SNPs.sig)
head(sig_SNPs.sig)

#################################
#New cols for reporting matches
#################################

#new col in gff table for snp counts
nrow(gff)
gff$sig_SNP_count = rep(0, nrow(gff))
#new cols in SNP table for whether a gene was found for a snp and type of relationship (no.match, in.gene, nearest), gene id of match, and distance of match (if nearest neighbor)
sig_SNPs.sig$match.type = vector(mode = "character", length = nrow(sig_SNPs.sig))
sig_SNPs.sig$geneID = vector(mode = "character", length = nrow(sig_SNPs.sig))
sig_SNPs.sig$distance = vector(mode = "numeric", length = nrow(sig_SNPs.sig))

#set to data.table to using rolling count update
gff.dt = as.data.table(gff)

#################################################

#Loop through SNPs
#Need to add gene name matches to SNP table (also deal with case of >1 match)
#Also add matched SNP count to genes 
for(i in 1:nrow(sig_SNPs.sig)){
    
    #filter for a gene that contains the SNP within start and stop positions
    temp = gff.dt[
        contig == sig_SNPs.sig$scaffold[i] &
        start <=  sig_SNPs.sig$position[i] &
        stop >= sig_SNPs.sig$position[i]
    ]
    #check if there is a matching gene. If not go to next
    if(nrow(temp) == 1){ #need the if in cases where there is no in-gene match
        gff.dt[ #the first lines are dt conditionals
            contig == sig_SNPs.sig$scaffold[i] &
                start <=  sig_SNPs.sig$position[i] &
                stop >= sig_SNPs.sig$position[i],
            #rolling count update
            sig_SNP_count := sig_SNP_count + 1
        ] 
    
        sig_SNPs.sig$match.type[i] = "in.gene"
        sig_SNPs.sig$geneID[i] = temp$geneID
        sig_SNPs.sig$distance[i] = 0
    }else{ #if there is no within gene match filter to the nearest gene(s)
        #filter to the correct scaffold
        temp2 = gff.dt[
            contig == sig_SNPs.sig$scaffold[i] 
        ] 
        #calculate distance to genes and return the nearest neighbor(s)
        temp2.match = nearest_gene(sig_SNPs.sig$position[i], temp2)
        
        #Need to account for possible ties in distance. First look for single match
        if(nrow(temp2.match) == 1){
            print(nrow(temp2.match))
            gff.dt[
                geneID == temp2.match$geneID, #conditional filter to gene ID
                sig_SNP_count := sig_SNP_count + 1 #update count
            ]  
      
            sig_SNPs.sig$match.type[i] = "nearest"
            sig_SNPs.sig$geneID[i] = temp2.match$geneID
            sig_SNPs.sig$distance[i] = temp2.match$dist
        }else if(nrow(temp2.match) > 0){ #if there are no matches we will update below
            print(nrow(temp2.match))
            gff.dt[
                geneID == temp2.match$geneID, #conditional filter to gene ID
                sig_SNP_count := sig_SNP_count + 1 #update count
            ]  
            
            sig_SNPs.sig$match.type[i] = "nearest.multiple"
            sig_SNPs.sig$geneID[i] = temp2.match$geneID
            sig_SNPs.sig$distance[i] = temp2.match$dist
        }else{
            sig_SNPs.sig$match.type[i] = "no.match"
            sig_SNPs.sig$geneID[i] = "no.match"
            sig_SNPs.sig$distance[i] = 0
        }
    }#end alternative to within gene
}#end for loop
            
sig_SNPs.sig$match.type %>% unique
sig_SNPs.sig$distance %>% abs %>% mean #4818.864
sig_SNPs.sig$distance %>% abs %>% median #4117.5
sig_SNPs.sig %>% filter(distance > 0 | distance < 0) %>% pull(distance) %>% abs %>% range #56 12871
sig_SNPs.sig %>% filter(match.type == "in.gene") %>% nrow #3

gff.dt %>% filter(sig_SNP_count > 0) #19 associated with SNPs
gff.dt %>% filter(sig_SNP_count > 1) #7 have >1 SNP

#######################
#plots of snp distance
#######################
library(ggplot2)
p1 = ggplot(sig_SNPs.sig, aes(x = distance)) +
    geom_histogram(breaks = seq(-25000, 25000, 250)) +
    theme_bw() +
    labs(x = "Gene distance from SNP (bp)", y = "Number SNPs (bin width 250 bp)")
p1

pdf(paste0("figures/GxE/RDA/",set_me,".SNP_dist.pdf"))
p1
dev.off()

##############################
#Write results tables
##############################
write.table(gff.dt %>% filter(sig_SNP_count > 0), paste0("data/Nf_LFMM_tables/", set_me , ".gene_SNP_hits.txt"), col.names = T, row.names = F, sep = "\t", quote = F)
write.table(sig_SNPs.sig, "data/Nf_LFMM_tables/tmin.SNPs.gene_found.nearest_neighbors.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(gff.dt %>% filter(sig_SNP_count > 0) %>% select(geneID), "data/Nf_LFMM_tables/tmin.geneIDs.nearest_neighbors.txt", col.names = F, row.names = F, sep = "\t", quote = F)

