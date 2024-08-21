library(dplyr)
library(ggplot2)
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
coords.Nf$lat[coords.Nf$collection_period == "early"] = NA
coords.Nf$lon[coords.Nf$collection_period == "early"] = NA

coords.Nd$lat[coords.Nd$collection_period == "early"] = NA
coords.Nd$lon[coords.Nd$collection_period == "early"] = NA

#calcualte geographic distance
Nf.Dgeo <- dist(dismo::Mercator(coords.Nf[,c("lon", "lat")]))
Nd.Dgeo <- dist(dismo::Mercator(coords.Nd[,c("lon", "lat")]))
Nc.Dgeo <- dist(dismo::Mercator(coords.Nc[,c("lon", "lat")]))

#########################
#long format for plotting
Nf.Dgeo.long = reshape2::melt(Nf.Dgeo %>% as.matrix)
Nf.Dgeo.long %>% filter(Var1 == "NG2" & Var2 == "NG1") %>% nrow
Nd.Dgeo.long = reshape2::melt(Nd.Dgeo %>% as.matrix)
Nc.Dgeo.long = reshape2::melt(Nc.Dgeo %>% as.matrix)
#Need to set self comps to NA
Nf.Dgeo.long[Nf.Dgeo.long$Var1 == Nf.Dgeo.long$Var2, "value"] = NA
Nd.Dgeo.long[Nd.Dgeo.long$Var1 == Nd.Dgeo.long$Var2, "value"] = NA
Nc.Dgeo.long[Nc.Dgeo.long$Var1 == Nc.Dgeo.long$Var2, "value"] = NA

#set gen dists to NA where geo dist == NA
Nf.Dgen.long = reshape2::melt(dist.Nf %>% as.matrix)
Nd.Dgen.long = reshape2::melt(dist.Nd %>% as.matrix)
Nc.Dgen.long = reshape2::melt(dist.Nc %>% as.matrix)

Nf.Dgen.long[is.na(Nf.Dgeo.long$value), "value"] = NA
Nd.Dgen.long[is.na(Nd.Dgeo.long$value), "value"] = NA
Nc.Dgen.long[is.na(Nc.Dgeo.long$value), "value"] = NA

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
Nf.Dgen.long$difsPerKb = Nf.Dgen.long$value / (41018940 / 1000)
Nd.Dgen.long$difsPerKb = Nd.Dgen.long$value / (38535154 / 1000)
Nc.Dgen.long$difsPerKb = Nc.Dgen.long$value / (40630626 / 1000)

#counts of instances
sum(Nf.Dgen.long$comp == "within", na.rm = T)
#682
sum(Nf.Dgen.long$comp == "between", na.rm = T)
#11974
sum(Nd.Dgen.long$comp == "within", na.rm = T)
#40
sum(Nd.Dgen.long$comp == "between", na.rm = T)
#560
Nc.Dgen.long[Nc.Dgen.long$comp == "within",]
sum(Nc.Dgen.long$comp == "within", na.rm = T)
#6
sum(Nc.Dgen.long$comp == "between", na.rm = T)
#14

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

#Nf and Nd shared within site comps
Nf_Nd_shared_sites = unique(Nd.withins.site$site)[unique(Nd.withins.site$site) %in% unique(Nf.withins.site$site)]

comp.withins = rbind(Nf.withins.site, Nd.withins.site, Nc.withins.site)
nrow(comp.withins)
#between
Nf.betweens = Nf.Dgen.long[Nf.Dgen.long$comp == "between" & !is.na(Nf.Dgen.long$value),]
Nf.betweens$spp = "Nf"
Nd.betweens = Nd.Dgen.long[Nd.Dgen.long$comp == "between" & !is.na(Nd.Dgen.long$value),]
Nd.betweens$spp = "Nd"
Nc.betweens = Nc.Dgen.long[Nc.Dgen.long$comp == "between" & !is.na(Nc.Dgen.long$value),]
Nc.betweens$spp = "Nc"
comp.betweens = rbind(Nf.betweens, Nd.betweens, Nc.betweens)

nrow(comp.betweens)
Nf_Nd_site_within = comp.withins %>% filter(site %in%  Nf_Nd_shared_sites)


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
#plot
Nf.Dgen.long$spp = "Nf"
Nd.Dgen.long$spp = "Nd"
Nc.Dgen.long$spp = "Nc"
all.Dgen = rbind(Nf.Dgen.long, Nd.Dgen.long, Nc.Dgen.long)
all.Dgen = all.Dgen[!is.na(all.Dgen$value),]

within_between_spp.aov = aov(difsPerKb ~ comp * spp, data = all.Dgen)
qqnorm(residuals(within_between_spp.aov))
plot(residuals(within_between_spp.aov))
summary(within_between_spp.aov)
TukeyHSD(within_between_spp.aov)


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
n_tab = data.frame(
    lab = paste("n sites =", c(Nf.within.site_n, Nd.within.site_n, Nc.within.site_n)),
    spp = c("Nf", "Nd", "Nc"),
    y = c(50, 7, 4)
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
    scale_x_continuous(limits = c(0,15)) +
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
    geom_text(data = n_tab, x = 12, aes(y = y, label = lab)) +
    my_gg_theme +
    labs(x = "Hamming distance (SNPs per Kb)", y = "Within site comparisons (count)", title = "a") +
    theme(
        plot.title = element_text(hjust = -0.12)
    )
p1
#just between
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
    scale_x_continuous(limits = c(0,15)) +
    my_gg_theme +
    labs(x = "Hamming distance (SNPs per Kb)", y = "PBetween site comparisons (count)", title = "b") +
    theme(
        plot.title = element_text(hjust = -0.12)
    )

pdf("figures/pop_gen/IBD/within-between.pdf", width = 16, height = 7)
grid.arrange(p1,p2, ncol = 2)
dev.off()
