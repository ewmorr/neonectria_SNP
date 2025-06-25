library(vegan)

pixy.Nf = readRDS("data/Nf/pixy/contig_all_dif_pop/dxy_dist.rds") 
pixy.Nf
class(pixy.Nf)

pixy.Nd = readRDS("data/Nd/pixy/contig_all_dif_pop/dxy_dist.rds") 
head(pixy.Nd)

pixy.Nc = readRDS("data/Nc/pixy/contig_all_dif_pop/dxy_dist.rds") 
head(pixy.Nc)

#get order of samples
dist.order.Nf = as.matrix(pixy.Nf ) %>% rownames
dist.order.Nd = as.matrix(pixy.Nd ) %>% rownames 
dist.order.Nc = as.matrix(pixy.Nc ) %>% rownames


dist.Nf = read.table("data/Nf/final_tables/rm_dups/FINAL_invariant.IBD_analyses.dist", header = F) 
dist.Nf = dist.Nf/41040857
dist.ID.Nf = read.table("data/Nf/final_tables/rm_dups/FINAL_invariant.IBD_analyses.dist.id", header = F)
rownames(dist.Nf) = dist.ID.Nf[,1]
colnames(dist.Nf) = dist.ID.Nf[,2]
#order
dist.Nf = dist.Nf[dist.order.Nf, dist.order.Nf]
dist.Nf = dist.Nf %>% as.dist()
dist.Nf

dist.Nd = read.table("data/Nd/final_tables/rm_dups/FINAL_invariant.IBD_analyses.dist", header = F) 
dist.Nd = dist.Nd/38535154
dist.ID.Nd = read.table("data/Nd/final_tables/rm_dups/FINAL_invariant.IBD_analyses.dist.id", header = F)
rownames(dist.Nd) = dist.ID.Nd[,1]
colnames(dist.Nd) = dist.ID.Nd[,2]
#order
dist.Nd = dist.Nd[dist.order.Nd, dist.order.Nd]
dist.Nd = dist.Nd %>% as.dist()
head(dist.Nd)

dist.Nc = read.table("data/Nc/final_tables/FINAL_invariant.IBD_analyses.dist", header = F) 
dist.Nc = dist.Nc/40630626
dist.ID.Nc = read.table("data/Nc/final_tables/FINAL_invariant.IBD_analyses.dist.id", header = F)
rownames(dist.Nc) = dist.ID.Nc[,1]
colnames(dist.Nc) = dist.ID.Nc[,2]
#order
dist.Nc = dist.Nc[dist.order.Nc, dist.order.Nc]
dist.Nc = dist.Nc %>% as.dist()
head(dist.Nc)

#mantel test
mantel(dist.Nf, pixy.Nf)
#Mantel statistic r: 0.9994 
#      Significance: 0.001 
mantel(dist.Nd, pixy.Nd)
#Mantel statistic r: 0.9979 
#      Significance: 0.001 
mantel(dist.Nc, pixy.Nc)
#Mantel statistic r:     1 
#      Significance: 0.0083333 

#long format for plotting
Nf.pixy.long = reshape2::melt(pixy.Nf %>% as.matrix)
Nd.pixy.long = reshape2::melt(pixy.Nd %>% as.matrix)
Nc.pixy.long = reshape2::melt(pixy.Nc %>% as.matrix)

Nf.dist.long = reshape2::melt(dist.Nf %>% as.matrix)
Nd.dist.long = reshape2::melt(dist.Nd %>% as.matrix)
Nc.dist.long = reshape2::melt(dist.Nc %>% as.matrix)


plot( Nf.pixy.long$value ~ Nf.dist.long$value)
plot( Nd.pixy.long$value ~ Nd.dist.long$value)
plot( Nc.pixy.long$value ~ Nc.dist.long$value)
