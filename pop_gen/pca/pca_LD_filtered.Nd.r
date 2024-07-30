library(pcadapt)
library(ggplot2)
library(dplyr)
source("R_scripts/ggplot_theme.txt")
library(RColorBrewer)
library(gridExtra)
source("R_scripts/make_site_metadata.r")

makeBed <- "data/Nd_SPANDx_all_seqs/out.filtered.LD_filtered_0.5_10Kb.bed"
makeVcf = "data/Nd_SPANDx_all_seqs/out.filtered.LD_filtered_0.5_10Kb.vcf"
makeBed.file <- read.pcadapt(makeBed, type = "bed")
bed_mat.makeBed = pcadapt::bed2matrix(makeBed.file)

sum(bed_mat.makeBed[is.na(bed_mat.makeBed) == F])

sum(bed_mat.makeBed[is.na(bed_mat.makeBed) == F] == 0)
sum(bed_mat.makeBed[is.na(bed_mat.makeBed) == F] == 1)
sum(bed_mat.makeBed[is.na(bed_mat.makeBed) == F] == 2)

x <- pcadapt(input = makeBed.file, K = 15, ploidy = 2) #K must be less than the number of species
#Ploidy is set to 2 because the function throws an error otherwise. This looks like because the conversion from BED is not recognizing haploid data. A look at `bed_mat.makeBed` shows SNPs encoded as 0 or 2 with no outcomes of 1 (i.e., het gentypes).


plot(x, option = "screeplot")

plot(x, option = "scores")
plot(x, option = "manhattan") #This is useful for full dataset
plot(x, option = "qqplot") #Another representation of p-value for SNP
plot(x, option = "stat.distribution") #chi.sq distribution

#Extract data for pretty plotting

#Plot scree plot (can use any of the different K val PCAs
#The singular values vector within PCA object contains the square root of porportion variance explained
scree_dat = data.frame(prop.var = x$singular.values^2, PC = 1:length(x$singular.values))

p1 = ggplot(scree_dat, aes(PC, prop.var)) +
geom_point() +
geom_line() +
scale_x_continuous(breaks = c(1,10,20)) +
labs(x = "PC axis", y = "Proportion variance explained") +
my_gg_theme
p1

pdf("figures/Nd.LD_filtered.PCA_scree_plot.pdf", width = 6, height = 5)
p1
dev.off()

#Plot PCA scores
sample_metadata
fam_info


pc_scores.metadata = data.frame(
    PC1 = x$scores[,1],
    PC2 = x$scores[,2],
    PC3 = x$scores[,3],
    PC4 = x$scores[,4],
    PC5 = x$scores[,5],
    PC6 = x$scores[,6],
    PC7 = x$scores[,7],
    PC8 = x$scores[,8],
    sample = fam_info.Nd[,1]
)

pc_scores.metadata = left_join(pc_scores.metadata, sample_metadata, by = "sample")

#25 color palette
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
pie(rep(1, 25), col = c25)

p1 = ggplot(pc_scores.metadata, aes(PC1, PC2, color = state.name)) +
geom_point(size = 3) +
scale_color_manual(values = c25) +
guides(color = "none") +
labs(
    x =
    paste(
            paste("PC1 (",round(scree_dat[1,1]*100, 1), sep = ""
            ),
        "% variance)", sep = ""
    ),
    y =
    paste(
        paste("PC2 (",round(scree_dat[2,1]*100, 1), sep = ""
        ),
        "% variance)", sep = ""
    ),
) +
my_gg_theme
p1

p1.guide = ggplot(pc_scores.metadata, aes(PC1, PC2, color = state.name)) +
geom_point(size = 3) +
scale_color_manual(values = c25) +
labs(
    x =
    paste(
            paste("PC1 (",round(scree_dat[1,1]*100, 1), sep = ""
            ),
        "% variance)", sep = ""
    ),
    y =
    paste(
        paste("PC2 (",round(scree_dat[2,1]*100, 1), sep = ""
        ),
        "% variance)", sep = ""
    ),
) +
my_gg_theme
p1.guide

pdf("figures/Nd.LD_filtered.PCA_plot.P1-P2.pdf", width = 8, height = 5)
p1.guide
dev.off()

p2 = ggplot(pc_scores.metadata, aes(PC1, PC3, color = state.name)) +
geom_point(size = 3) +
scale_color_manual(values = c25) +
guides(color = "none") +
labs(
    x =
    paste(
            paste("PC1 (",round(scree_dat[1,1]*100, 1), sep = ""
            ),
        "% variance)", sep = ""
    ),
    y =
    paste(
        paste("PC3 (",round(scree_dat[3,1]*100, 1), sep = ""
        ),
        "% variance)", sep = ""
    ),
) +
my_gg_theme
p2

p3 = ggplot(pc_scores.metadata, aes(PC2, PC3, color = state.name)) +
geom_point(size = 3) +
scale_color_manual(values = c25) +
labs(
    x =
    paste(
            paste("PC2 (",round(scree_dat[2,1]*100, 1), sep = ""
            ),
        "% variance)", sep = ""
    ),
    y =
    paste(
        paste("PC3 (",round(scree_dat[3,1]*100, 1), sep = ""
        ),
        "% variance)", sep = ""
    ),
) +
my_gg_theme 
p3


pdf("figures/Nd.LD_filtered.PCA_plot.pdf", width = 16, height = 4)
grid.arrange(p1,p2,p3, nrow = 1, widths = c(0.3, 0.3, 0.4))
#grid.arrange(p1,p2,p3,p4,p5,p6, nrow = 2, widths = c(0.3, 0.3, 0.4))
dev.off()

#For full dataset will want to rerun PCA, determine number of axes of pop structure, and then recalculate with that number K. Then test statistics for outlier SNPs can be assessed

c25 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "palegreen2", "#FDBF6F", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4", "brown")
pie(rep(1, 25), col = c25)
