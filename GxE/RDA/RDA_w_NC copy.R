library(vegan)
library(dplyr)
library(ggplot2)
library(GGally)
source("library/ggplot_theme.txt")


Y = read.table("data/Nf/GxE/RDA/Y.wNC.filteredMAC.tsv", header = F, sep = "\t")
SNP_pos = read.csv("data/Nf/GxE/RDA/SNP_pos.wNC.filteredMAC.csv")
sample_ids = read.table("data/Nf/GxE/RDA/sample_ids.wNC", header = F)

#read site metadata
site.info = read.csv("data/sample_metadata/site_info.csv")
site.bioclim = read.csv("data/sample_metadata/Nf.site_bioclim.csv")
site.bioclim = site.bioclim %>% filter(state != "WV")
site.bioclim
source("data/sample_metadata/world_clim/bioclim_names.txt")
bioclim_var_names

############################
# setting up metadata
site_metadata = left_join(
    site.bioclim,
    site.info %>% select(state, lat, lon, duration_infection),
    by = c("state", "lat", "lon")
) %>% unique() 

#sample metadata. Needs to be sorted by sampleIDs 
sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
rownames(sample_metadata.Nf) = sample_metadata.Nf$Sequence_label

sample_metadata.Nf.sorted = sample_metadata.Nf[sample_ids$V1,]
nrow(sample_metadata.Nf.sorted)
sample_metadata.site_info = left_join(sample_metadata.Nf.sorted, site_metadata)
is.na(sample_metadata.site_info$bio1) %>% sum()

row_ids = which(sample_ids$V1 %in% sample_metadata.site_info$Sequence_label)

head(sample_metadata.site_info)

X = sample_metadata.site_info %>%
    select(starts_with("bio"), duration_infection)
X

###########################################
# examining which vars need to be excluded
# 
#scale the vars
X.uniq = X %>% unique()
nrow(X.uniq)

site_metadata.scaled = apply(
    X.uniq, 
    2, 
    scale
) %>% as.data.frame

varcors = cor(site_metadata.scaled)
sort(rowSums(abs(varcors) > 0.7))
bioclim_var_names
# we can start with removing bio1 (MAT), warmest/coldest quarter (bio10, bio11), MAP (bio12) first
# we keep bio6 for now (min temp coldest month) bc it may be biologically important
site_metadata.scaled.reduced = site_metadata.scaled %>%
    select(-bio1, -bio10, -bio11, -bio12)
varcors = cor(site_metadata.scaled.reduced)
sort(rowSums(abs(varcors) > 0.7))
bioclim_var_names
# lets remove the precip quarter estimates (bio16,17)
# these are repettive with the -est month measures
# remove bio16, bio17
site_metadata.scaled.reduced = site_metadata.scaled.reduced %>%
    select(-bio16, -bio17)
varcors = cor(site_metadata.scaled.reduced)
sort(rowSums(abs(varcors) > 0.7))
bioclim_var_names[names(rowSums(varcors > 0.7))]
# remove temp seasonality (bio4), bio9

site_metadata.scaled.reduced = site_metadata.scaled.reduced %>%
    select(-bio4, -bio9) 
varcors = cor(site_metadata.scaled.reduced)
sort(rowSums(abs(varcors) > 0.7))
bioclim_var_names[names(rowSums(varcors > 0.7))]

site_metadata.scaled.reduced = site_metadata.scaled.reduced %>%
    select(-bio7) 
varcors = cor(site_metadata.scaled.reduced)
sort(rowSums(abs(varcors) > 0.7))

site_metadata.scaled.reduced = site_metadata.scaled.reduced %>%
    select(-bio19) 
varcors = cor(site_metadata.scaled.reduced)
sort(rowSums(abs(varcors) > 0.7))

site_metadata.scaled.reduced = site_metadata.scaled.reduced %>%
    select(-bio13, -bio14) 
varcors = cor(site_metadata.scaled.reduced)
sort(rowSums(abs(varcors) > 0.7))
bioclim_var_names[names(rowSums(varcors > 0.7))]



p1 = ggpairs(
    site_metadata.scaled.reduced, 
    columns = 1:ncol(site_metadata.scaled.reduced)
) +
    my_gg_theme.def_size +
    theme(
        axis.text.x = element_text(angle = 55, hjust = 1)
    )
p1

bioclim_var_names[names(rowSums(varcors > 0.7))]


pdf("figures/GxE/RDA/env_vars/pairs_cor.w_NC.pdf", width = 16, height = 16)
p1
dev.off()

names(varcors)

X[,rownames(varcors)]

# let's try using a pre-se;ected set (not based on cors)
selected_env = c("bio2", "bio4", "bio5", "bio6", "bio8", "bio9", "bio13", "bio14", "bio18", "bio19", "duration_infection")

#https://popgen.nescent.org/2018-03-27_RDA_GEA.html
#
###################
# there is something of an art to picking the X matrix still
# 
# reduced vars
env_mat = X[,rownames(varcors)]
# full vars
env_mat = X
#Some constraints or conditions were aliased because they were redundant. This can happen if terms are linearly dependent (collinear): ‘bio7’

# selected vars
env_mat = X[,selected_env]

############
# run the RDA
Y.rda <- rda(Y~ ., data=env_mat, scale=T)
Y.rda

#saveRDS(Y.rda, "data/Nf/GxE/RDA/w_NC.RDA_env_selected.rds")

RsquareAdj(Y.rda)

#The eigenvalues for the constrained axes reflect the variance explained by each canonical axis:
summary(eigenvals(Y.rda, model = "constrained"))
# We can visualize this information using a screeplot of the canonical eigenvalues by calling screeplot:
screeplot(Y.rda)

# scree
ggplot() +
geom_line(aes(x=c(1:length(Y.rda$CCA$eig)), y=as.vector(Y.rda$CCA$eig)/sum(Y.rda$CCA$eig)*100), linetype="dotted",
size = 1.5, color="darkgrey") +
geom_point(aes(x=c(1:length(Y.rda$CCA$eig)), y=as.vector(Y.rda$CCA$eig)/sum(Y.rda$CCA$eig)*100), size = 3,
color="darkgrey") +
scale_x_continuous(name = "Ordination axes", breaks=c(1:10)) +
ylab("% variance") +
theme_bw()

# reduced = 3 (but not a strong decline in var explained)
# selected = 3 (could make an arg for 4, but it's quite a bit less)
# full 5

#Now let’s check our RDA model for significance using formal tests. We can assess both the full model and each constrained axis using F-statistics (Legendre et al, 2010). 
signif.full <- anova.cca(Y.rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full

#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: rda(formula = Y ~ bio2 + bio3 + bio5 + bio6 + bio8 + bio15 + bio18 + duration_infection, data = env_mat, scale = T)
#         Df Variance      F Pr(>F)    
#Model     8    11673 1.2735  0.001 ***
#Residual 89   101978                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#########################
# The full model is significant, but that doesn’t tell us much. We can check each constrained axis for significance using the code below. For this test, each constrained axis is tested using all previous constrained axes as conditions. See ?anova.cca and Legendre et al. (2010) for details. The purpose here is to determine which constrained axes we should investigate for candidate loci.

# This analysis is time intensive (taking up to a few hours for the full wolf data set), so we will not run the code here. If we did run it, we would find that the first three constrained axes are significant (p = 0.001); constrained axis 4 has a p-value of 0.080, while axes 5-8 have p-values > 0.850. This corresponds with our evaluation of the screeplot, above.

###
### This analysis is also mem intensive. Pulls in about 31 G on the server

#start.time <- Sys.time()
#signif.axis <- anova.cca(Y.rda, by="axis", parallel=getOption("mc.cores"))
#end.time <- Sys.time()
#time.taken <- end.time - start.time
#time.taken

signif.axis = readRDS("data/Nf/GxE/RDA/w_NC.signif_axis.env_reduced.rds")
signif.axis
# 3 sig
signif.axis = readRDS("data/Nf/GxE/RDA/w_NC.signif_axis.env_selected.rds")
signif.axis
#  sig
signif.axis = readRDS("data/Nf/GxE/RDA/w_NC.signif_axis.env_full.rds")
signif.axis
#  sig

#########################


# vegan has a simple function for checking Variance Inflation Factors for the predictor variables used in the model:
vif.cca(Y.rda)
#              bio2               bio3               bio5               bio8               bio9              bio13              bio14              bio18    duration_infection
#          9.699805          16.352897           7.642433           8.283968          11.560292          11.859774           5.657037          10.049121    3.366768   
          

#All values are below 10, and most are below 5, which indicates that multicollinearity among these predictors shouldn’t be a problem for the model. We could remove one of the temperature variables (AMT or MDR) if we were concerned about these higher VIF values (Zuur et al., 2010).

# 
# 
#              bio1               bio2               bio3               bio4               bio5               bio6               bio7               bio8               bio9 
#      6.257379e+05       1.211844e+04       3.245226e+04       3.178912e+05       2.041892e+04       5.478910e+04                 NA       2.736696e+03       5.257181e+01 
#             bio10              bio11              bio12              bio13              bio14              bio15              bio16              bio17              bio18 
#      4.583412e+04       1.819367e+06       2.713571e+04       4.614652e+03       1.589045e+04       9.553935e+03       6.755058e+02       1.160403e+04       7.318879e+03 
#             bio19 duration_infection 
#      1.151517e+03       5.033007e+01 

#selected vars
#      bio2       bio4       bio5       bio6       bio8       bio9      bio13      bio14      bio18      bio19 
# 21.194469 156.824176  59.167557 320.658452   9.044443  13.622518   9.102534  24.790246  10.329052  26.959872 

#######################
#######################
#
# We’ll start with simple triplots from vegan. Here we’ll use scaling=3 (also known as “symmetrical scaling”) for the ordination plots. This scales the SNP and individual scores by the square root of the eigenvalues. See Borcard et al. (2011) or the vegan help for more information on scaling in RDA plots.

plot(Y.rda, scaling=3)
plot(Y.rda, choices = c(1, 3), scaling=3)
plot(Y.rda, choices = c(1, 4), scaling=3)
plot(Y.rda, choices = c(1, 5), scaling=3)
bioclim_var_names

#We’ll use the loadings of the SNPs in the ordination space to determine which SNPs are candidates for local adaptation. The SNP loadings are stored as species in the RDA object. We’ll extract the SNP loadings from the three significant constrained axes:
load.rda <- scores(Y.rda, choices=c(1:4), display="species")  # Species scores for the first three constrained axes
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA1")
hist(load.rda[,4], main="Loadings on RDA2")

# get outliers based on two-tailed z score
# from Forester
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}


nrow(load.rda)

#selected env
cand1 <- outliers(load.rda[,1],3) # 787
cand2 <- outliers(load.rda[,2],3) # 287
cand3 <- outliers(load.rda[,3],3) # 746
cand4 <- outliers(load.rda[,4],3) # 387
length(cand1)
length(cand2)
length(cand3)
length(cand4)


(length(cand1) + length(cand2)) - length( intersect(names(cand1), names(cand2)) ) #can not use setdiff bc it only returns first vec
length( intersect(names(cand1), names(cand2)) )
#0 are shared
length( intersect(names(cand1), names(cand3)) )
# 1
length( intersect(names(cand1), names(cand4)) )
# 1
length( intersect(names(cand2), names(cand3)) )
# 0
length( intersect(names(cand2), names(cand4)) )
# 0
length( intersect(names(cand3), names(cand4)) )
# 5


#from Capblanq, following pcadapt strategy
rdadapt<-function(rda,K)
{
    zscores<-Y.rda$CCA$v[,1:as.numeric(K)]
    resscale <- apply(zscores, 2, scale)
    resmaha <- robust::covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
    lambda <- median(resmaha)/qchisq(0.5,df=K)
    reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
    qval <- qvalue::qvalue(reschi2test)
    q.values_rdadapt<-qval$qvalues
    return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}
# env selected
# 3 axes reduced
# 3 axes selected 
# 3 axes full (could make an arg for 5)
last_axis = 5

res_rdadapt<-rdadapt(Y.rda, last_axis)

#note that capblanq also uses q < 0.1
which(res_rdadapt$q.values < 0.05) %>% length
# full data set 5 axes had 3131 at 0.05
# selected vars has 1403 at 0.05
# non-correalted vars 3 axes has 1355
which(res_rdadapt$q.values < 0.1) %>% length
# selected vars including duration infection 2039 at 0.1 with 3 axes
# full vars 5 axes 4551
# non-correalted vars 3 axes has 1853
 
# manhattan plot outliers
qvalue_max = 0.1

ggplot() +
geom_point(
    aes(
        x=c(1:length(res_rdadapt[,1])), 
        y=-log10(res_rdadapt[,1])
    ),
    col = "gray83"
) +
geom_point(
    aes(
        x=c(1:length(res_rdadapt[,1]))[which(res_rdadapt[,2] <qvalue_max)],
        y = -log10(res_rdadapt[which(res_rdadapt[,2] < qvalue_max),1])
    ),
    col = "orange"
) +
xlab("SNPs") + ylab("-log10(p.values)") +
theme_bw()

####################################
# projection of outliers in rda space
# function for axis choice
axis_choice_biplot = function(rda.obj = NULL, choices=c(1,2), res_rdadapt = res_rdadapt, qvalue_max = qvalue_max){
    rda.scores = scores(rda.obj, choices = choices, scaling = "symmetric")
    
    ggplot() +
        geom_point(
            aes(
                x=rda.scores$species[which(res_rdadapt$q.values >= qvalue_max),1], 
                y=rda.scores$species[which(res_rdadapt$q.values >= qvalue_max),2]
            ), 
            col = "gray86"
        ) +
        geom_point(
            aes(
                x=rda.scores$species[which(res_rdadapt$q.values < qvalue_max),1], 
                y=rda.scores$species[which(res_rdadapt$q.values < qvalue_max),2]
            ), 
            col = "orange"
        ) +
        geom_segment(
            aes(
                xend=rda.scores$biplot[,1], 
                yend=rda.scores$biplot[,2], 
                x=0, 
                y=0
            ),
            colour="black", size=0.5, linetype=1,
            arrow=arrow(length = unit(0.02, "npc"))
        ) +
        geom_text(
            aes(
                x=1.125*rda.scores$biplot[,1], 
                y=1.125*rda.scores$biplot[,2], 
                label = rownames(rda.scores$biplot)
            )
        ) +
    xlab(paste("RDA",choices[1])) + ylab(paste("RDA",choices[2])) +
    theme_bw() +
    theme(legend.position="none")
}
###################################################################

axis_choice_biplot(rda.obj = Y.rda, choices = c(1,2), res_rdadapt = res_rdadapt, qvalue_max = qvalue_max)
axis_choice_biplot(rda.obj = Y.rda, choices = c(1,3), res_rdadapt = res_rdadapt, qvalue_max = qvalue_max)
axis_choice_biplot(rda.obj = Y.rda, choices = c(1,4), res_rdadapt = res_rdadapt, qvalue_max = qvalue_max)
axis_choice_biplot(rda.obj = Y.rda, choices = c(1,5), res_rdadapt = res_rdadapt, qvalue_max = qvalue_max)

axis_choice_biplot(rda.obj = Y.rda, choices = c(2,3), res_rdadapt = res_rdadapt, qvalue_max = qvalue_max)
axis_choice_biplot(rda.obj = Y.rda, choices = c(2,4), res_rdadapt = res_rdadapt, qvalue_max = qvalue_max)
axis_choice_biplot(rda.obj = Y.rda, choices = c(2,5), res_rdadapt = res_rdadapt, qvalue_max = qvalue_max)

# plotting genome pos
nrow(SNP_pos) == nrow(res_rdadapt)
SNP_pos$rdadapt_qvalue = res_rdadapt$q.values
SNP_pos$rdadapt_sig = ifelse(SNP_pos$rdadapt_qvalue < qvalue_max, "sig", "n.s.")

ggplot(
    SNP_pos %>% filter(length > 100000), 
    aes(x = position/10^6, y = rdadapt_qvalue, color = rdadapt_sig)
) +
#facet_wrap(~scaffold) +
facet_grid(. ~ scaffold, scales = "free_x", space='free') +
geom_point(alpha = 0.5, size = 1) +
scale_x_continuous(breaks = c(0, seq(from = 1, to = 6, by = 1)) ) +
scale_y_continuous(
    trans  = scales::compose_trans("log10", "reverse"),
    labels = scales::label_log()

) +
scale_color_manual(values = c("grey", "black"), guide = "none") +
my_gg_theme +
labs(x = "Position (Mbp)", y = "") +
theme(
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
)

##################################
# using the correlation based approach from forester
# with the SNPs identified using the prdadapt method

dim(Y)
Y.rdadapt_sig = Y[,which(res_rdadapt$q.values < qvalue_max)]
SNP_pos.rdadapt_sig = SNP_pos[which(res_rdadapt$q.values < qvalue_max),]

# cors with environment

snp_env_cors = matrix(nrow = ncol(Y.rdadapt_sig), ncol =  ncol(env_mat))
colnames(snp_env_cors) = colnames(env_mat)

for (i in 1:ncol(Y.rdadapt_sig)) {
  snp_env_cors[i,] <- apply(env_mat,2,function(x) cor(x,Y.rdadapt_sig[,i]))
}

SNP_pos.rdadapt_sig = cbind(SNP_pos.rdadapt_sig, snp_env_cors)
head(SNP_pos.rdadapt_sig)

SNP_pos.rdadapt_sig$max_env_cor = vector(mode = "numeric", length = nrow(SNP_pos.rdadapt_sig))
SNP_pos.rdadapt_sig$max_env_predictor = vector(mode = "character", length = nrow(SNP_pos.rdadapt_sig))

for (i in 1:nrow(SNP_pos.rdadapt_sig)) {
  temp.cors <- snp_env_cors[i,]
  SNP_pos.rdadapt_sig$max_env_predictor[i] = names(which.max(abs(temp.cors))) # gives the variable
  SNP_pos.rdadapt_sig$max_env_cor[i] = temp.cors[which.max(abs(temp.cors))]       # gives the correlation
}

last_col = ncol(SNP_pos.rdadapt_sig)
SNP_pos.rdadapt_sig[,c(1,2,last_col-1,last_col)]

# cors with RDA axes
#last_axis = 5
rda.axes = data.frame(scores(Y.rda, choices = 1:last_axis, display = "sites"))
snp_axes_cors = matrix(nrow = ncol(Y.rdadapt_sig), ncol =  ncol(rda.axes))
colnames(snp_axes_cors) = colnames(rda.axes)

for (i in 1:ncol(Y.rdadapt_sig)) {
  snp_axes_cors[i,] <- apply(rda.axes,2,function(x) cor(x,Y.rdadapt_sig[,i]))
}

SNP_pos.rdadapt_sig = cbind(SNP_pos.rdadapt_sig, snp_axes_cors)
SNP_pos.rdadapt_sig$max_rdaAxis_cor = vector(mode = "numeric", length = nrow(SNP_pos.rdadapt_sig))
SNP_pos.rdadapt_sig$max_rdaAxis_predictor = vector(mode = "character", length = nrow(SNP_pos.rdadapt_sig))

for (i in 1:nrow(SNP_pos.rdadapt_sig)) {
  temp.cors <- snp_axes_cors[i,]
  SNP_pos.rdadapt_sig$max_rdaAxis_predictor[i] = names(which.max(abs(temp.cors))) # gives the variable
  SNP_pos.rdadapt_sig$max_rdaAxis_cor[i] = temp.cors[which.max(abs(temp.cors))]       # gives the correlation
}

new_last_col = ncol(SNP_pos.rdadapt_sig)
SNP_pos.rdadapt_sig[,c(1,2,3,last_col-1,last_col,new_last_col-1,new_last_col)]

bioclim_var_names

table(SNP_pos.rdadapt_sig$max_env_predictor)
table(SNP_pos.rdadapt_sig$max_rdaAxis_predictor)

# selected vars
#             bio13              bio14              bio18              bio19               bio5               bio6               bio9 
#                 5                  6                  8                  5                  6                  3                  2 
# duration_infection 
#                 1 

#reduced vars
#             bio14              bio18               bio5               bio8 duration_infection 
#                61                  1                  1                  2                  1 

#full vars
#              bio1              bio10              bio11              bio12              bio13              bio14 
#               151                457                 14                275                429                377 
#             bio15              bio16              bio17              bio18              bio19               bio2 
#                96                237                 45                165                 19                 88 
#              bio3               bio4               bio5               bio6               bio7               bio8 
#                41                 20                122                 43                131                115 
#              bio9 duration_infection 
#                26                242 
                
#write.csv(SNP_pos.rdadapt_sig, "data/Nf/GxE/RDA/w_NC.rdadapt_full_env.csv", quote = F, row.names = F)
#write.csv(SNP_pos.rdadapt_sig, "data/Nf/GxE/RDA/w_NC.rdadapt_selected_env.csv", quote = F, row.names = F)
#write.csv(SNP_pos.rdadapt_sig, "data/Nf/GxE/RDA/w_NC.rdadapt_reduced_env.csv", quote = F, row.names = F)
    
colnames(SNP_pos.rdadapt_sig)
    SNP_pos.rdadapt_sig[,c(1,2,3,19,20)]
    
    