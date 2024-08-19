library(dplyr)
library(ggplot2)
source("library/ggplot_theme.txt")

dist.raw.Nf = readRDS("data/Nf/IBD/hamming_dist.rds")
dist.raw.Nd = readRDS("data/Nd/IBD/hamming_dist.rds")
dist.raw.Nc = readRDS("data/Nc/IBD/hamming_dist.rds")
#get order of samples
dist.order.Nf = as.matrix(dist.raw.Nf ) %>% rownames 
dist.order.Nd = as.matrix(dist.raw.Nd ) %>% rownames 
dist.order.Nc = as.matrix(dist.raw.Nc ) %>% rownames 

length(dist.order.Nf)
length(dist.order.Nd)
length(dist.order.Nc)

#metadata
sample_metadata.Nf = read.csv("data/sample_metadata/Nf_filtered.lat_lon_dur_inf.csv")
sample_metadata.Nd = read.csv("data/sample_metadata/Nd_filtered.lat_lon_dur_inf.csv")
nrow(sample_metadata.Nf)
nrow(sample_metadata.Nd)
rownames(sample_metadata.Nf) = sample_metadata.Nf$Sequence_label
rownames(sample_metadata.Nd) = sample_metadata.Nd$Sequence_label
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

#For now we only want to compare isolates collected in the modern era
# We may use the older isolates for a different test
coords.Nf$lat[coords.Nf$collection_period == "early"] = NA
coords.Nf$lon[coords.Nf$collection_period == "early"] = NA

coords.Nd$lat[coords.Nd$collection_period == "early"] = NA
coords.Nd$lon[coords.Nd$collection_period == "early"] = NA

#calcualte geographic distance
Nf.Dgeo <- dist(dismo::Mercator(coords.Nf[,c("lon", "lat")]))
Nd.Dgeo <- dist(dismo::Mercator(coords.Nd[,c("lon", "lat")]))

#########################
#long format for plotting
Nf.Dgeo.long = reshape2::melt(Nf.Dgeo %>% as.matrix)
Nd.Dgeo.long = reshape2::melt(Nd.Dgeo %>% as.matrix)
#Need to set self comps to NA
Nf.Dgeo.long[Nf.Dgeo.long$Var1 == Nf.Dgeo.long$Var2, "value"] = NA
Nd.Dgeo.long[Nd.Dgeo.long$Var1 == Nd.Dgeo.long$Var2, "value"] = NA

#set gen dists to NA where geo dist == NA
Nf.Dgen.long = reshape2::melt(dist.raw.Nf %>% as.matrix)
Nd.Dgen.long = reshape2::melt(dist.raw.Nd %>% as.matrix)

Nf.Dgen.long[is.na(Nf.Dgeo.long$value), "value"] = NA
Nd.Dgen.long[is.na(Nd.Dgeo.long$value), "value"] = NA

#create within between cats
Nf.Dgen.long$comp = ifelse(Nf.Dgeo.long$value == 0, "within", "between")
Nd.Dgen.long$comp = ifelse(Nd.Dgeo.long$value == 0, "within", "between")

#convert distances to nucelotide difs (instead of proportional) and
# correct for genome size (see IBD.individuals.calculate_hamming for number comps)
Nf.Dgen.long$difsPerKb = Nf.Dgen.long$value * 999132 #/ 42948.211 
Nd.Dgen.long$difsPerKb = Nd.Dgen.long$value * 1413869 #/ 44950.817
range(Nf.Dgen.long$difsPerKb, na.rm = T)
range(Nd.Dgen.long$difsPerKb, na.rm = T)

#counts of instances
sum(Nf.Dgen.long$comp == "within", na.rm = T)
#682
sum(Nf.Dgen.long$comp == "between", na.rm = T)
#11974
sum(Nd.Dgen.long$comp == "within", na.rm = T)
#40
sum(Nd.Dgen.long$comp == "between", na.rm = T)
#560

#dfs for anova
#within
Nf.withins = Nf.Dgen.long[Nf.Dgen.long$comp == "within" & !is.na(Nf.Dgen.long$value),]
Nf.withins$spp = "Nf"
Nd.withins = Nd.Dgen.long[Nd.Dgen.long$comp == "within" & !is.na(Nd.Dgen.long$value),]
Nd.withins$spp = "Nd"
comp.withins = rbind(Nf.withins, Nd.withins)
nrow(comp.withins)
#between
Nf.betweens = Nf.Dgen.long[Nf.Dgen.long$comp == "between" & !is.na(Nf.Dgen.long$value),]
Nf.betweens$spp = "Nf"
Nd.betweens = Nd.Dgen.long[Nd.Dgen.long$comp == "between" & !is.na(Nd.Dgen.long$value),]
Nd.betweens$spp = "Nd"
comp.betweens = rbind(Nf.betweens, Nd.betweens)
nrow(comp.betweens)

#aov within sites
aov.withins = aov(difsPerKb ~ spp, comp.withins)
qqnorm(residuals(aov.withins))
plot(residuals(aov.withins))
summary(aov.withins)
#aov between sites
aov.betweens = aov(difsPerKb ~ spp, comp.betweens)
qqnorm(residuals(aov.betweens))
plot(residuals(aov.betweens))
summary(aov.betweens)

#plot
Nf.Dgen.long$spp = "Nf"
Nd.Dgen.long$spp = "Nd"
all.Dgen = rbind(Nf.Dgen.long, Nd.Dgen.long)
all.Dgen = all.Dgen[!is.na(all.Dgen$value),]

ggplot(all.Dgen,
       aes(x = difsPerKb)
) +
    geom_histogram() +
    facet_grid(spp~comp, scales = "free_y") +
    my_gg_theme

#just wihtin
ggplot(all.Dgen %>% filter(comp == "within"),
       aes(x = difsPerKb)
) +
    geom_histogram(binwidth = 0.1) +
    facet_wrap(
        ~factor(spp, levels = c("Nf", "Nd"), 
                labels = c("N. faginata", "N. ditissima")
        ),
        scales = "free_y", 
        ncol = 1
    ) +
    my_gg_theme +
    labs(x = "Hamming distance")
#just between
ggplot(all.Dgen %>% filter(comp == "between"),
       aes(x = difsPerKb)
) +
    geom_histogram(binwidth = 0.1) +
    facet_wrap(
        ~factor(spp, levels = c("Nf", "Nd"), 
                labels = c("N. faginata", "N. ditissima")
        ),
        scales = "free_y", 
        ncol = 1
    ) +
    my_gg_theme +
    labs(x = "Hamming distance")

