library(dplyr)
library(ggplot2)
library(geosphere)
library(gridExtra)
source("library/ggplot_theme.txt")

dist.Nf = read.table("data/Nf/final_tables/rm_dups/FINAL_invariant.IBD_analyses.dist", header = F) 
dist.ID.Nf = read.table("data/Nf/final_tables/rm_dups/FINAL_invariant.IBD_analyses.dist.id", header = F)
rownames(dist.Nf) = dist.ID.Nf[,1]
colnames(dist.Nf) = dist.ID.Nf[,2]
dist.Nf = dist.Nf %>% as.dist()
dist.Nf

dist.Nd = read.table("data/Nd/final_tables/rm_dups/FINAL_invariant.IBD_analyses.dist", header = F) 
dist.ID.Nd = read.table("data/Nd/final_tables/rm_dups/FINAL_invariant.IBD_analyses.dist.id", header = F)
rownames(dist.Nd) = dist.ID.Nd[,1]
colnames(dist.Nd) = dist.ID.Nd[,2]
dist.Nd = dist.Nd %>% as.dist()
head(dist.Nd)

dist.Nc = read.table("data/Nc/final_tables/FINAL_invariant.IBD_analyses.dist", header = F) 
dist.ID.Nc = read.table("data/Nc/final_tables/FINAL_invariant.IBD_analyses.dist.id", header = F)
rownames(dist.Nc) = dist.ID.Nc[,1]
colnames(dist.Nc) = dist.ID.Nc[,2]
dist.Nc = dist.Nc %>% as.dist()
head(dist.Nc)

#get order of samples
dist.order.Nf = as.matrix(dist.Nf ) %>% rownames
dist.order.Nd = as.matrix(dist.Nd ) %>% rownames 
dist.order.Nc = as.matrix(dist.Nc ) %>% rownames

length(dist.order.Nf)
length(dist.order.Nd)
length(dist.order.Nc)

#metadata
sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
sample_metadata.Nd = read.csv("data/sample_metadata/Nd_filtered.lat_lon_dur_inf.csv")
sample_metadata.Nc = read.csv("data/sample_metadata/Nc_canton_loc_date.lat_lon.csv")

nrow(sample_metadata.Nf)
nrow(sample_metadata.Nd)
nrow(sample_metadata.Nc)
rownames(sample_metadata.Nf) = sample_metadata.Nf$Sequence_label
rownames(sample_metadata.Nd) = sample_metadata.Nd$Sequence_label
rownames(sample_metadata.Nc) = sample_metadata.Nc$Sequence_label


#lat lon tables (ordered)
coords.Nf = data.frame(
    samp = sample_metadata.Nf[dist.order.Nf, "Sequence_label"],
    collection_period = sample_metadata.Nf[dist.order.Nf, "collection_period"],
    lat = sample_metadata.Nf[dist.order.Nf, "lat"],
    lon = sample_metadata.Nf[dist.order.Nf, "lon"],
    stringsAsFactors = T
)
rownames(coords.Nf) = coords.Nf$samp

coords.Nd = data.frame(
    samp = sample_metadata.Nd[dist.order.Nd, "Sequence_label"],
    collection_period = sample_metadata.Nd[dist.order.Nd, "collection_period"],
    lat = sample_metadata.Nd[dist.order.Nd, "lat"],
    lon = sample_metadata.Nd[dist.order.Nd, "lon"],
    stringsAsFactors = T
)
rownames(coords.Nd) = coords.Nd$samp

coords.Nc = data.frame(
    samp = sample_metadata.Nc[dist.order.Nc, "Sequence_label"],
    collection_period = rep("modern", nrow(sample_metadata.Nc)),
    lat = sample_metadata.Nc[dist.order.Nc, "lat"],
    lon = sample_metadata.Nc[dist.order.Nc, "lon"],
    stringsAsFactors = T
)
rownames(coords.Nc) = coords.Nc$samp

#For now we only want to compare isolates collected in the modern era
# We may use the older isolates for a different test

#coords.Nf$lat[coords.Nf$collection_period == "early"] = NA
#coords.Nf$lon[coords.Nf$collection_period == "early"] = NA

#coords.Nd$lat[coords.Nd$collection_period == "early"] = NA
#coords.Nd$lon[coords.Nd$collection_period == "early"] = NA

#calcualte geographic distance
#old bad way #dist(dismo::Mercator(coords.Nf[,c("lon", "lat")]))
Nf.Dgeo <- distm(x = coords.Nf[,c("lon", "lat")], fun = distVincentyEllipsoid)
Nd.Dgeo <- distm(x = coords.Nd[,c("lon", "lat")], fun = distVincentyEllipsoid)
Nc.Dgeo <- distm(x = coords.Nc[,c("lon", "lat")], fun = distVincentyEllipsoid)

#########################
#long format for plotting
#we need to first NA the upper triangle of the matrix to avoid repeat comps
#just do this here and use these to set NAs in the gen dist mat
Nf.mat = as.matrix(Nf.Dgeo)
Nf.mat[upper.tri(Nf.mat, diag=T)] = NA
Nf.Dgeo.long = reshape2::melt(Nf.mat)
Nf.Dgeo.long %>% filter(Var1 == "NG2" & Var2 == "NG1") %>% nrow
Nd.mat = as.matrix(Nd.Dgeo)
Nd.mat[upper.tri(Nd.mat, diag=T)] = NA
Nd.Dgeo.long = reshape2::melt(Nd.mat)
Nc.mat = as.matrix(Nc.Dgeo)
Nc.mat[upper.tri(Nc.mat, diag=T)] = NA
Nc.Dgeo.long = reshape2::melt(Nc.mat)
#Need to set self comps to NA
##don't need this anymore bc did it in the matrix convert above
#Nf.Dgeo.long[Nf.Dgeo.long$Var1 == Nf.Dgeo.long$Var2, "value"] = NA
#Nd.Dgeo.long[Nd.Dgeo.long$Var1 == Nd.Dgeo.long$Var2, "value"] = NA
#Nc.Dgeo.long[Nc.Dgeo.long$Var1 == Nc.Dgeo.long$Var2, "value"] = NA

#set gen dists to NA where geo dist == NA
Nf.Dgen.long = reshape2::melt(dist.Nf %>% as.matrix)
Nd.Dgen.long = reshape2::melt(dist.Nd %>% as.matrix)
Nc.Dgen.long = reshape2::melt(dist.Nc %>% as.matrix)

Nf.Dgen.long[is.na(Nf.Dgeo.long$value), "value"] = NA
Nd.Dgen.long[is.na(Nd.Dgeo.long$value), "value"] = NA
Nc.Dgen.long[is.na(Nc.Dgeo.long$value), "value"] = NA

#now remove NA values from everything
Nf.Dgeo.long = Nf.Dgeo.long[!is.na(Nf.Dgeo.long$value),]
Nd.Dgeo.long = Nd.Dgeo.long[!is.na(Nd.Dgeo.long$value),]
Nc.Dgeo.long = Nc.Dgeo.long[!is.na(Nc.Dgeo.long$value),]
Nf.Dgen.long = Nf.Dgen.long[!is.na(Nf.Dgen.long$value),]
Nd.Dgen.long = Nd.Dgen.long[!is.na(Nd.Dgen.long$value),]
Nc.Dgen.long = Nc.Dgen.long[!is.na(Nc.Dgen.long$value),]

nrow(Nf.Dgeo.long)
nrow(Nf.Dgen.long)

#create within between cats
Nf.Dgen.long$comp = ifelse(Nf.Dgeo.long$value == 0, "within", "between")
Nd.Dgen.long$comp = ifelse(Nd.Dgeo.long$value == 0, "within", "between")
Nc.Dgen.long$comp = ifelse(Nc.Dgeo.long$value == 0, "within", "between")

#convert raw nucelotide difs from plink to per Kb corrected for genome size 
# i.e., / the number of sites in each table / 1000
# #see the filtration.sh scripts for each species for number of sites
# 
# Nf L741 FINAL_invariant.nuclear.vcf.gz = 41018940
# Nd L391 FINAL_invariant.nuclear.vcf.gz = 38535154
# Nc L274 FINAL_invariant.nuclear.vcf.gz = 40630626 bp
Nf.Dgen.long$difsPerKb = Nf.Dgen.long$value / (41040857 / 1000)
Nd.Dgen.long$difsPerKb = Nd.Dgen.long$value / (38535154 / 1000)
Nc.Dgen.long$difsPerKb = Nc.Dgen.long$value / (40630626 / 1000)

#counts of instances
sum(Nf.Dgen.long$comp == "within", na.rm = T)
#341
sum(Nf.Dgen.long$comp == "between", na.rm = T)
#5987
#6214 with early isolates
sum(Nd.Dgen.long$comp == "within", na.rm = T)
#20
#23 with early
sum(Nd.Dgen.long$comp == "between", na.rm = T)
#280
#412 with early
Nc.Dgen.long[Nc.Dgen.long$comp == "within",]
sum(Nc.Dgen.long$comp == "within", na.rm = T)
#3
sum(Nc.Dgen.long$comp == "between", na.rm = T)
#7

#dfs for anova
#within
Nf.withins = Nf.Dgen.long[Nf.Dgen.long$comp == "within" & !is.na(Nf.Dgen.long$value),]
Nf.withins$spp = "Nf"
Nf.join_site = data.frame(
    Var1 = sample_metadata.Nf$Sequence_label, site = sample_metadata.Nf$state
)
Nf.withins.site = left_join(Nf.withins, Nf.join_site)
Nf.within.site_n = Nf.withins.site$site %>% unique %>% length
Nd.withins = Nd.Dgen.long[Nd.Dgen.long$comp == "within" & !is.na(Nd.Dgen.long$value),]
Nd.withins$spp = "Nd"
Nd.join_site = data.frame(
    Var1 = sample_metadata.Nd$Sequence_label, site = sample_metadata.Nd$state
)
Nd.withins.site = left_join(Nd.withins, Nd.join_site)
Nd.within.site_n = Nd.withins.site$site %>% unique %>% length

Nc.withins = Nc.Dgen.long[Nc.Dgen.long$comp == "within" & !is.na(Nc.Dgen.long$value),]
Nc.withins$spp = "Nc"
Nc.join_site = data.frame(
    Var1 = sample_metadata.Nc$Sequence_label, site = sample_metadata.Nc$Canton
)
Nc.withins.site = left_join(Nc.withins, Nc.join_site)
Nc.within.site_n = Nc.withins.site$site %>% unique %>% length

#bind the species tables
comp.withins = rbind(Nf.withins.site, Nd.withins.site, Nc.withins.site)
nrow(comp.withins)
#364
#367 with early

#between
Nf.betweens = Nf.Dgen.long[Nf.Dgen.long$comp == "between" & !is.na(Nf.Dgen.long$value),]
Nf.betweens$spp = "Nf"
Nf.join_site1 = data.frame(
    Var1 = sample_metadata.Nf$Sequence_label, site1 = sample_metadata.Nf$state
)
Nf.join_site2 = data.frame(
    Var2 = sample_metadata.Nf$Sequence_label, site2 = sample_metadata.Nf$state
)
Nf.betweens.site = left_join(Nf.betweens, Nf.join_site1) %>%
    left_join(., Nf.join_site2)
Nf.betweens.site.n = paste(Nf.betweens.site$site1, Nf.betweens.site$site2, sep = "-") %>% 
    unique() %>% 
    length()

Nd.betweens = Nd.Dgen.long[Nd.Dgen.long$comp == "between" & !is.na(Nd.Dgen.long$value),]
Nd.betweens$spp = "Nd"
Nd.join_site1 = data.frame(
    Var1 = sample_metadata.Nd$Sequence_label, site1 = sample_metadata.Nd$state
)
Nd.join_site2 = data.frame(
    Var2 = sample_metadata.Nd$Sequence_label, site2 = sample_metadata.Nd$state
)
Nd.betweens.site = left_join(Nd.betweens, Nd.join_site1) %>%
    left_join(., Nd.join_site2)
Nd.betweens.site.n = paste(Nd.betweens.site$site1, Nd.betweens.site$site2, sep = "-") %>% 
    unique() %>% 
    length()

Nc.betweens = Nc.Dgen.long[Nc.Dgen.long$comp == "between" & !is.na(Nc.Dgen.long$value),]
Nc.betweens$spp = "Nc"
Nc.join_site1 = data.frame(
    Var1 = sample_metadata.Nc$Sequence_label, site1 = sample_metadata.Nc$Canton
)
Nc.join_site2 = data.frame(
    Var2 = sample_metadata.Nc$Sequence_label, site2 = sample_metadata.Nc$Canton
)
Nc.betweens.site = left_join(Nc.betweens, Nc.join_site1) %>%
    left_join(., Nc.join_site2)
Nc.betweens.site.n = paste(Nc.betweens.site$site1, Nc.betweens.site$site2, sep = "-") %>% 
    unique() %>% 
    length()

#bind spp tables
comp.betweens = rbind(Nf.betweens, Nd.betweens, Nc.betweens)
nrow(comp.betweens)
#6274
# 6633 with early

#################
#aov within sites
aov.withins = aov(difsPerKb ~ spp, comp.withins)
qqnorm(residuals(aov.withins))
plot(residuals(aov.withins))
summary(aov.withins)
TukeyHSD(aov.withins)

#lme within sites
library(nlme)
comp.withins$spp = comp.withins$spp %>% as.factor
comp.withins$site = comp.withins$site %>% as.factor
lme.withins = lme(fixed = difsPerKb ~ spp, random = ~1|site, data = comp.withins)
qqnorm(residuals(lme.withins))
plot(residuals(lme.withins))
summary(lme.withins)
library(multcomp)
summary(glht(lme.withins, linfct = mcp(spp = "Tukey")))

#aov between sites
aov.betweens = aov(difsPerKb ~ spp, comp.betweens)
qqnorm(residuals(aov.betweens))
plot(residuals(aov.betweens))
summary(aov.betweens)
TukeyHSD(aov.betweens)

#single table with all species and both comps
Nf.Dgen.long$spp = "Nf"
Nd.Dgen.long$spp = "Nd"
Nc.Dgen.long$spp = "Nc"
all.Dgen = rbind(Nf.Dgen.long, Nd.Dgen.long, Nc.Dgen.long)
all.Dgen = all.Dgen[!is.na(all.Dgen$value),]

#two-way ANOVA comparing within/between and spp
within_between_spp.aov = aov(difsPerKb ~ comp * spp, data = all.Dgen)
qqnorm(residuals(within_between_spp.aov))
plot(residuals(within_between_spp.aov))
summary(within_between_spp.aov)
TukeyHSD(within_between_spp.aov)


#compare Nf and Nd only within the same site
#Nf and Nd shared within site comps
Nf_Nd_shared_sites = unique(Nd.withins.site$site)[unique(Nd.withins.site$site) %in% unique(Nf.withins.site$site)]
Nf_Nd_site_within = comp.withins %>% filter(site %in%  Nf_Nd_shared_sites)
Nf_Nd_site_within$site %>% unique()
Nf_Nd_site_within.aov = aov(difsPerKb ~ site/spp, Nf_Nd_site_within)
Nf_Nd_site_within.aov
summary(Nf_Nd_site_within.aov)
TukeyHSD(Nf_Nd_site_within.aov)

#planned comparisons (i.e., spp within site)
Nf_Nd_site_within.aov = aov(difsPerKb ~ spp, Nf_Nd_site_within %>% filter(site == "WV"))
qqnorm(residuals(Nf_Nd_site_within.aov))
plot(residuals(Nf_Nd_site_within.aov))
summary(Nf_Nd_site_within.aov)
#0.000161 ***
Nf_Nd_site_within.aov = aov(difsPerKb ~ spp, Nf_Nd_site_within %>% filter(site == "MI"))
qqnorm(residuals(Nf_Nd_site_within.aov))
plot(residuals(Nf_Nd_site_within.aov))
summary(Nf_Nd_site_within.aov)
#not enough data (only one obs per)
Nf_Nd_site_within.aov = aov(difsPerKb ~ spp, Nf_Nd_site_within %>% filter(site == "MI.UP"))
qqnorm(residuals(Nf_Nd_site_within.aov))
plot(residuals(Nf_Nd_site_within.aov))
summary(Nf_Nd_site_within.aov)
#not enough data (only one obs per)
Nf_Nd_site_within.aov = aov(difsPerKb ~ spp, Nf_Nd_site_within %>% filter(site == "ME.N"))
qqnorm(residuals(Nf_Nd_site_within.aov))
plot(residuals(Nf_Nd_site_within.aov))
summary(Nf_Nd_site_within.aov)
#1.68e-07 ***
#sig where possible to compare
################################
################################
#Bootstrap 95% CIs
#
head(all.Dgen)

spps = all.Dgen$spp %>% unique()
comps = all.Dgen$comp %>% unique()
options(warn=2)

boot_df = data.frame(
    spp = rep(c("Nf", "Nd", "Nc"), 2),
    comp = c(rep("within",3), rep("between", 3)),
    med.boot = vector(mode = "numeric", length = 6),
    lower.CI = vector(mode = "numeric", length = 6),
    upper.CI = vector(mode = "numeric", length = 6)
)

for(i in spps){
    for(j in comps){
        temp.Dgen = all.Dgen %>% filter(spp == i & comp == j)
        temp.Dgen[!is.na(temp.Dgen$difsPerKb),]
        temp.boot = replicate(
            9999, 
            mean(
                sample(
                    temp.Dgen[!is.na(temp.Dgen$difsPerKb),"difsPerKb"], 
                    replace = T
                )
            )
        )
        quants.CI = quantile(temp.boot, probs = c(0.025, 0.5, 0.975))
        boot_df[boot_df$spp == i & boot_df$comp == j, "med.boot"] = quants.CI[2]
        boot_df[boot_df$spp == i & boot_df$comp == j, "lower.CI"] = quants.CI[1]
        boot_df[boot_df$spp == i & boot_df$comp == j, "upper.CI"] = quants.CI[3]
    }
}
options(warn=1)

boot_df
# without early
# 
#  spp    comp   med.boot  lower.CI   upper.CI
#1  Nf  within  3.5619717  3.470350  3.6439682
#2  Nd  within  9.2735259  9.056061  9.4319463
#3  Nc  within  0.7778024  0.768002  0.7964116
#4  Nf between  4.0208136  4.012703  4.0287666
#5  Nd between  9.2552116  9.168552  9.3158062
#6  Nc between 13.3665252 12.637448 13.8768721

# with early
# 
#  spp    comp   med.boot  lower.CI   upper.CI
#1  Nf  within  3.5612472  3.470086  3.6411180
#2  Nd  within  9.2234957  9.024972  9.3784264
#3  Nc  within  0.7778024  0.768002  0.7964116
#4  Nf between  4.0248340  4.016795  4.0328803
#5  Nd between  9.1931483  9.132322  9.2393101
#6  Nc between 13.3657376 12.642678 13.8749629

################
################
#plots
#
ggplot(all.Dgen,
       aes(x = difsPerKb)
) +
    geom_histogram(binwidth = 0.1) +
    facet_grid(spp~comp, scales = "free_y") +
    my_gg_theme

#################
#################
# annotate each panel with the number of sites (within) or site pairs (bn)
# Third panel with bootstrapped means and 95% CIs
# 
#just within
#
n_tab.within = data.frame(
    lab = paste("n sites =", c(Nf.within.site_n, Nd.within.site_n, Nc.within.site_n)),
    spp = c("Nf", "Nd", "Nc"),
    y = c(70, 9, 5)
)
p1 = ggplot(all.Dgen %>% filter(comp == "within"),
       aes(x = difsPerKb)
) +
    geom_histogram(binwidth = 0.1) +
    facet_wrap(
        ~factor(spp, levels = c("Nf", "Nd", "Nc"), 
                labels = c("N. faginata", "N. ditissima", "N. coccinea")
        ),
        scales = "free_y", 
        ncol = 1
    ) +
    coord_cartesian(xlim = c(0,15)) +
    #ggh4x::scale_y_facet(
    #    spp == "Nd", scale_y_continuous(breaks = c(0,5,10))
    #) +
    ggh4x::facetted_pos_scales(
        y = list(
            #spp == "Nf" ~ scale_y_continuous(limits = c(0, 80)),
            #"spp" == "Nd" ~ scale_y_continuous(breaks = c(0,5,10))
            #spp == "Nc" ~ scale_y_continuous(trans = "reverse")
            scale_y_continuous(breaks = c(0, 40, 80)),
            scale_y_continuous(breaks = c(0, 5, 10)),
            scale_y_continuous(breaks = c(0, 3, 6))
        )
    ) +
    geom_text(data = n_tab.within, x = 12.5, size = 5, hjust = 0, aes(y = y, label = lab)) +
    my_gg_theme +
    labs(x = "Hamming distance (SNPs per Kb)", y = "Comparisons of individuals\nwithin sites (count)", title = "a") +
    theme(
        plot.title = element_text(hjust = -0.12)
    )
p1
#just between
n_tab.between = data.frame(
lab = paste("n site pairs =", c(Nf.betweens.site.n, Nd.betweens.site.n, Nc.betweens.site.n)),
spp = c("Nf", "Nd", "Nc"),
y = c(1420, 55, 6.7)
)

p2 = ggplot(all.Dgen %>% filter(comp == "between"),
       aes(x = difsPerKb)
) +
    geom_histogram(binwidth = 0.1) +
    facet_wrap(
        ~factor(spp, levels = c("Nf", "Nd", "Nc"), 
                labels = c("N. faginata", "N. ditissima", "N. coccinea")
        ),
        scales = "free_y", 
        ncol = 1
    ) +
    coord_cartesian(xlim = c(0,15)) +
    geom_text(data = n_tab.between, x = -0.5, size = 5, aes(y = y, label = lab), hjust = 0) +
    my_gg_theme +
    labs(x = "Hamming distance (SNPs per Kb)", y = "Comparisons of individuals\nbetween sites (count)", title = "b") +
    theme(
        plot.title = element_text(hjust = -0.12)
    )
p2

pdf("figures/pop_gen/IBD/within-between.include_early.pdf", width = 16, height = 7)
grid.arrange(p1,p2, ncol = 2)
dev.off()

#plot all grid with ggh4x::facet_grid2

n_tab.wnbn = data.frame(
    lab = c(
        paste("sites n =", c(Nf.within.site_n, Nd.within.site_n, Nc.within.site_n)),
        paste("site pairs\nn =", c(Nf.betweens.site.n, Nd.betweens.site.n, Nc.betweens.site.n))
    ),
    spp = rep(c("Nf", "Nd", "Nc"), 2),
    y = c(40.5, 4.38, 2.625, 750, 40, 3.2),
    comp = c(rep("within", 3), rep("between", 3))
) #the y is nicely aligned between cols at 7 inch height pdf

p1 = ggplot(all.Dgen, 
       aes(x = difsPerKb)
) +
    geom_histogram(binwidth = 0.1) +
    geom_text(data = n_tab.wnbn, x = c(rep(11, 3), rep(0, 3)), size = 3, aes(y = y, label = lab), hjust = 0) +
    ggh4x::facet_grid2(
        rows = vars(factor(
            spp, 
            levels = c("Nf", "Nd", "Nc")#, 
#            labels = c("N. faginata", "N. ditissima", "N. coccinea")
        )),
        cols = vars(factor(
            comp,
            levels = c("within", "between")#,
#            labels = c("p", "q") # for some reason this breaks when adding labels
        )),
        scales = "free_y", 
        independent = "y",
        switch = "y",
        labeller = labeller(
            .rows = sppNames <- c("Nf" = "N. faginata", "Nd" = "N. ditissima", "Nc" = "N. coccinea"),
            .cols = compNames <- c("within" = "Within site comparisons", "between" = "Between site comparisons")
        )
    ) +    
    ggh4x::facetted_pos_scales(
        y = list(
            scale_y_continuous(breaks = c(0, 15, 30, 45)),
            scale_y_continuous(breaks = c(0, 300, 600, 900)),
            scale_y_continuous(breaks = c(0, 2, 4, 6)),
            #scale_y_continuous(breaks = c(0, 10,20, 30)),
            scale_y_continuous(breaks = c(0, 15,30, 45)),
            scale_y_continuous(breaks = c(0, 1, 2, 3)),
            scale_y_continuous(breaks = c(0, 1,2,3,4))
        )
    ) +
    my_gg_theme.def_size +
    labs(
        x = expression(paste("Hamming distance (SNPs Kb"^-1,")")), 
        y = "Pairwise comparisons of individuals (count)", 
        title = "a"
    ) +
    theme(
        plot.title = element_text(hjust = -0.125, margin = margin(b = -7.5)),
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        strip.placement = "outside"
    )
p1
pdf("figures/pop_gen/IBD/within-between.grid.include_early.pdf", width = 7, height = 4)
p1
dev.off()


p2 = ggplot(boot_df, 
    aes(
        y = med.boot, 
        x = factor(spp,
            levels = c("Nc", "Nd", "Nf"), 
            labels = c("N. coccinea", "N. ditissima", "N. faginata")
        ),
        shape = factor(comp, levels = c("between", "within"))
    )
) +
    geom_point(position = position_dodge(width = 0.9), size = 2.5) +
    geom_errorbar(
        aes(ymin = lower.CI, ymax = upper.CI), 
        position = position_dodge(width = 0.9),
        width = 0.4
    ) +
    scale_shape_manual(values = c(1,2), breaks = c("within", "between")) +
    #geom_linerange(
    #    aes(ymin = lower.CI, ymax = upper.CI), 
    #    position = position_dodge(width = 0.9),
    #    linewidth = 3
    #) +
    my_gg_theme.def_size +
    coord_flip() +
    scale_y_continuous(breaks = c(0,5,10), limits = c(0,NA)) +
    labs(
        shape = "Site comparison", 
        y = expression(paste("Hamming distance (SNPs Kb"^-1,")")),
        title = "b"
    ) +
    theme(
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        legend.position = "top", #c(0.75,0.85),
        legend.margin = margin(0, 0, -10, 0),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = -0.075, margin = margin(b = -7.5))
    ) 
p2


pdf("figures/pop_gen/IBD/within-between.CIs.include_early.pdf", width = 10, height = 4)
#p1 + p2
grid.arrange(p1,p2,ncol = 2, widths = c(0.7,0.3))
dev.off()


p2 = ggplot(boot_df, 
            aes(
                y = med.boot, 
                x = factor(spp,
                           levels = c("Nc", "Nd", "Nf"), 
                           labels = c("N. coccinea", "N. ditissima", "N. faginata")
                ),
                shape = factor(comp, levels = c("between", "within")),
                color = factor(comp, levels = c("between", "within"))
            )
) +
    geom_point(position = position_dodge(width = 0.9), size = 3) +
    geom_errorbar(
        aes(ymin = lower.CI, ymax = upper.CI), 
        position = position_dodge(width = 0.9),
        width = 0.4,
        color = "black"
    ) +
    scale_shape_manual(values = c(16,17), breaks = c("within", "between")) +
    scale_color_brewer(palette = "Dark2", breaks = c("within", "between")) +
    #geom_linerange(
    #    aes(ymin = lower.CI, ymax = upper.CI), 
    #    position = position_dodge(width = 0.9),
    #    linewidth = 3
    #) +
    my_gg_theme.def_size +
    coord_flip() +
    scale_y_continuous(breaks = c(0,5,10), limits = c(0,NA)) +
    labs(
        shape = "Site comparison", 
        color = "Site comparison",
        y = expression(paste("Hamming distance (SNPs Kb"^-1,")")),
        title = "b"
    ) +
    theme(
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        legend.position = "top", #c(0.75,0.85),
        legend.margin = margin(0, 0, -10, 0),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = -0.075, margin = margin(b = -7.5))
    ) 
p2

pdf("figures/pop_gen/IBD/within-between.CIs.color.include_early.pdf", width = 10, height = 4)
#p1 + p2
grid.arrange(p1,p2,ncol = 2, widths = c(0.7,0.3))
dev.off()


##################
##################    
#3 panel abc fig
#
#pulling figs from above
#y = c(70, 9, 5, 1250, 49, 5.85),

n_tab.within = data.frame(
    lab = paste("sites n =", c(Nf.within.site_n, Nd.within.site_n, Nc.within.site_n)),
    spp = c("Nf", "Nd", "Nc"),
    y = c(70, 9, 5)
)
p1 = ggplot(all.Dgen %>% filter(comp == "within"),
            aes(x = difsPerKb)
) +
    geom_histogram(binwidth = 0.1) +
    facet_wrap(
        ~factor(spp, levels = c("Nf", "Nd", "Nc"), 
                labels = c("N. faginata", "N. ditissima", "N. coccinea")
        ),
        scales = "free_y", 
        ncol = 1
    ) +
    coord_cartesian(xlim = c(0,15)) +
    #ggh4x::scale_y_facet(
    #    spp == "Nd", scale_y_continuous(breaks = c(0,5,10))
    #) +
    ggh4x::facetted_pos_scales(
        y = list(
            #spp == "Nf" ~ scale_y_continuous(limits = c(0, 80)),
            #"spp" == "Nd" ~ scale_y_continuous(breaks = c(0,5,10))
            #spp == "Nc" ~ scale_y_continuous(trans = "reverse")
            scale_y_continuous(breaks = c(0, 25, 50, 75)),
            scale_y_continuous(breaks = c(0, 4, 8)),
            scale_y_continuous(breaks = c(0, 2, 4, 6))
        )
    ) +
    geom_text(data = n_tab.within, x = 11, size = 3, hjust = 0, aes(y = y, label = lab)) +
    my_gg_theme.def_size +
    labs(x = "", y = "Comparisons of individuals\nwithin sites (count)", title = "a") +
    theme(
        plot.title = element_text(hjust = -0.25),
        strip.background = element_blank(),
        strip.text = element_blank()
    )
p1
#just between
n_tab.between = data.frame(
    lab = paste("site pairs\nn =", c(Nf.betweens.site.n, Nd.betweens.site.n, Nc.betweens.site.n)),
    spp = c("Nf", "Nd", "Nc"),
    y = c(1220, 49, 5.85)
)

p2 = ggplot(all.Dgen %>% filter(comp == "between"),
            aes(x = difsPerKb)
) +
    geom_histogram(binwidth = 0.1) +
    facet_wrap(
        ~factor(spp, levels = c("Nf", "Nd", "Nc"), 
                labels = c("N. faginata", "N. ditissima", "N. coccinea")
        ),
        scales = "free_y", 
        ncol = 1
    ) +
    coord_cartesian(xlim = c(0,15)) +
    geom_text(data = n_tab.between, x = 0, size = 3, aes(y = y, label = lab), hjust = 0) +
    my_gg_theme.def_size +
    labs(x = "Hamming distance (SNPs per Kb)", y = "between sites (count)", title = "b") +
    theme(
        plot.title = element_text(hjust = -0.25),
        strip.background = element_blank(),
        strip.text = element_blank()
    )
p2
p3 = ggplot(boot_df, 
            aes(
                x = med.boot, 
                y = factor(spp,
                           levels = c("Nc", "Nd", "Nf"), 
                           labels = c("N. coccinea", "N. ditissima", "N. faginata")
                ),
                shape = factor(comp, levels = c("within", "between"))
            )
) +
    geom_point(position = position_dodge(width = 0.9), size = 2.5) +
    geom_errorbar(
        aes(xmin = lower.CI, xmax = upper.CI), 
        position = position_dodge(width = 0.9),
        width = 0.4
    ) +
    facet_wrap(
        ~factor(spp, levels = c("Nf", "Nd", "Nc"), 
                labels = c("N. faginata", "N. ditissima", "N. coccinea")
        ),
        scales = "free_y", 
        strip.position = "right",
        ncol = 1
    ) +
    scale_shape_manual(values = c(1,2), labels = c("within" = "within sites", "between" = "between sites")) +
    coord_cartesian(xlim = c(0,15)) +
    #geom_linerange(
    #    aes(ymin = lower.CI, ymax = upper.CI), 
    #    position = position_dodge(width = 0.9),
    #    linewidth = 3
    #) +
    my_gg_theme.def_size +
    #coord_flip() +
    scale_x_continuous(breaks = c(0,5,10,15)) +
    labs(
        shape = "Site comparison", 
        y = "Bootstrap 95% CIs",
        x= "",
        title = "c"
    ) +
    theme(
        #axis.title.y = element_blank(),
        #axis.text.y = element_text(angle = 90, hjust = 0.5),
        #legend.position = c(0.75,0.85),
        legend.position = "bottom",
        legend.margin = margin(-26.5, 0, -2, 0),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = -0.12)
    ) 
p3

pdf("figures/pop_gen/IBD/within-between.three_panel.pdf", width = 10, height = 4)
grid.arrange(p1,p2,p3,ncol = 3)
dev.off()


