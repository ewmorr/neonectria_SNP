library(vegan)
library(dplyr)
library(ggplot2)
library(GGally)
source("library/ggplot_theme.txt")


Y = read.table("data/Nf/GxE/RDA/Y.noNC.filteredMAC.tsv", header = F, sep = "\t")
SNP_pos = read.csv("data/Nf/GxE/RDA/SNP_pos.noNC.filteredMAC.csv")
sample_ids = read.table("data/Nf/GxE/RDA/sample_ids.noNC", header = F)

#read site metadata
site.info = read.csv("data/sample_metadata/site_info.csv")
site.bioclim = read.csv("data/sample_metadata/Nf.site_bioclim.csv")
site.bioclim = site.bioclim %>% filter(state != "WV" & state != "NC")
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
# we can start with removing bio1, bio11, bio12 first
# we keep bio6 for now (min temp coldest month) bc it may be biologically important
site_metadata.scaled.reduced = site_metadata.scaled %>%
    select(-bio1, -bio11, -bio12)
varcors = cor(site_metadata.scaled.reduced)
sort(rowSums(abs(varcors) > 0.7))
bioclim_var_names
# lets remove the quarter estimates (wettest is one of the highest cors now)
# these are repettive with the -est month measures
# remove bio16, bio17, bio10
site_metadata.scaled.reduced = site_metadata.scaled.reduced %>%
    select(-bio16, -bio17, -bio10)
varcors = cor(site_metadata.scaled.reduced)
sort(rowSums(abs(varcors) > 0.7))
bioclim_var_names[names(rowSums(varcors > 0.7))]
# bio3 and 4 correlted with each other
# lets go with isothermaility (bio3)
# bio 5 with bio 6 only (0.7)
# bio 7 with bio 3,4,6,9
# let's go with -bio6
site_metadata.scaled.reduced = site_metadata.scaled.reduced %>%
    select(-bio6) 
varcors = cor(site_metadata.scaled.reduced)
sort(rowSums(abs(varcors) > 0.7))

site_metadata.scaled.reduced = site_metadata.scaled.reduced %>%
    select(-bio4, -bio7) 
varcors = cor(site_metadata.scaled.reduced)
sort(rowSums(abs(varcors) > 0.7))

site_metadata.scaled.reduced = site_metadata.scaled.reduced %>%
    select(-bio19) 
varcors = cor(site_metadata.scaled.reduced)
sort(rowSums(abs(varcors) > 0.7))

# last to deal with is precip seasonality versus precip direst month
# we also have precip wettest month included and these are easier to interpret 
# than how precip is distrbuted throughout the year
# we exclude bio 15
site_metadata.scaled.reduced = site_metadata.scaled.reduced %>%
    select(-bio15) 
varcors = cor(site_metadata.scaled.reduced)
sort(rowSums(abs(varcors) > 0.7))


p1 = ggpairs(
    site_metadata.scaled.reduced, 
    columns = 1:ncol(site_metadata.scaled.reduced)
) +
    my_gg_theme.def_size +
    theme(
        axis.text.x = element_text(angle = 55, hjust = 1)
    )

bioclim_var_names[names(rowSums(varcors > 0.7))]


pdf("figures/GxE/RDA/env_vars/pairs_cor.no_NC.pdf", width = 16, height = 16)
p1
dev.off()

names(varcors)

X[,rownames(varcors)]

# let's try using a pre-se;ected set (not based on cors)
selected_env = c("bio2", "bio4", "bio5", "bio6", "bio8", "bio9", "bio13", "bio14", "bio18", "bio19")

#https://popgen.nescent.org/2018-03-27_RDA_GEA.html
#
Y.rda <- rda(Y~ ., data=X[,rownames(varcors)], scale=T)
Y.rda <- rda(Y~ ., data=X, scale=T)
#Some constraints or conditions were aliased because they were redundant. This can happen if terms are linearly dependent (collinear): ‘bio7’
Y.rda <- rda(Y~ ., data=X[,selected_env], scale=T)
Y.rda

# saveRDS(Y.rda, "data/Nf/GxE/RDA/no_NC.RDA_env_selected.rds")

RsquareAdj(Y.rda)

#The eigenvalues for the constrained axes reflect the variance explained by each canonical axis:
summary(eigenvals(Y.rda, model = "constrained"))
# We can visualize this information using a screeplot of the canonical eigenvalues by calling screeplot:
screeplot(Y.rda)

#Now let’s check our RDA model for significance using formal tests. We can assess both the full model and each constrained axis using F-statistics (Legendre et al, 2010). 
signif.full <- anova.cca(Y.rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full

#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: rda(formula = Y ~ bio2 + bio3 + bio5 + bio8 + bio9 + bio13 + bio14 + bio18 + duration_infection, data = X[, rownames(varcors)], scale = T)
#         Df Variance     F Pr(>F)    
#Model     9    13959 1.215  0.001 ***
#Residual 82   104678                 
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

signif.axis = readRDS("data/Nf/GxE/RDA/no_NC.RDA.signif_axes.rds")
signif.axis

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



#######################
#######################
#plot
# We’ll start with simple triplots from vegan. Here we’ll use scaling=3 (also known as “symmetrical scaling”) for the ordination plots. This scales the SNP and individual scores by the square root of the eigenvalues. See Borcard et al. (2011) or the vegan help for more information on scaling in RDA plots.

plot(Y.rda, scaling=3)
plot(Y.rda, choices = c(1, 3), scaling=3)
plot(Y.rda, choices = c(1, 4), scaling=3)
bioclim_var_names




