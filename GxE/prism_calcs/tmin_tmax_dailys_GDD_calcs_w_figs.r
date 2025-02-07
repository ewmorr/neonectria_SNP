require(tidyverse)
require(RColorBrewer)
source("~/ggplot_theme.txt")

#saveRDS(sites_climate_list, "intermediate_RDS/sites_daily_tmin_tmax_20072018.list.rds")
sites_climate_list.df = readRDS(file = "sample_metadata/sites_daily_tmin_tmax_ppt_20072018.df.rds")

#dates
sites_climate_list.df$year = substr(sites_climate_list.df$day, 1, 4) %>% as.numeric
sites_climate_list.df$month = substr(sites_climate_list.df$day, 5, 6) %>% as.numeric
sites_climate_list.df$dom = substr(sites_climate_list.df$day, 7, 8) %>% as.numeric

sites_climate_list.df$date = ISOdate(
    sites_climate_list.df$year,
    sites_climate_list.df$month,
    sites_climate_list.df$dom,
    tz = "") %>%
as.Date

sites_climate_list.df$doy = strftime(sites_climate_list.df$date, format = "%j")

###############
#Calculate GDD#

#based 4 GDD (heating
sites_climate_list.df$HDD4 = (sites_climate_list.df$tmin + sites_climate_list.df$tmax)/2 - 4
#convert negatives to zero
sites_climate_list.df[, colnames(sites_climate_list.df) == "HDD4"][
    sites_climate_list.df[,colnames(sites_climate_list.df) == "HDD4"] < 0
] = 0

#base 20 CDD (chilling)
sites_climate_list.df$CDD20 = (273+(sites_climate_list.df$tmin + sites_climate_list.df$tmax)/2) - (273+20)
#convert positives to zero
sites_climate_list.df[,colnames(sites_climate_list.df) == "CDD20"][
    sites_climate_list.df[,colnames(sites_climate_list.df) == "CDD20"] > 0
] = 0

#convert all CDD before 8/15 to zero (based on Richardson et al. 2006)
#This does not work for leap year (Aug 14 still counted)
sites_climate_list.df[,colnames(sites_climate_list.df) == "CDD20"][
    sites_climate_list.df[,colnames(sites_climate_list.df) == "doy"] < 227
] = 0
#convert leap year days
sites_climate_list.df[,colnames(sites_climate_list.df) == "CDD20"][
    sites_climate_list.df[,colnames(sites_climate_list.df) == "month"] == 8 &
    sites_climate_list.df[,colnames(sites_climate_list.df) == "dom"] == 14
] = 0

#Also calculate daily freeze thaw events
sites_climate_list.df$freezeThaw = vector(mode = "numeric", length = length(sites_climate_list.df$tmin))

sites_climate_list.df[,colnames(sites_climate_list.df) == "freezeThaw"][
    sites_climate_list.df[,colnames(sites_climate_list.df) == "tmin"] <= 0 &
    sites_climate_list.df[,colnames(sites_climate_list.df) == "tmax"] > 0
] = 1

###################################
#Calculate cumulative sums of HDD4 and CDD20 in years/sites

sites_climate_list.df.cumsum = sites_climate_list.df %>% group_by(year, Site) %>% mutate(HDD4.cumsum = cumsum(HDD4))

sites_climate_list.df.cumsum = sites_climate_list.df.cumsum %>% group_by(year, Site) %>% mutate(CDD20.cumsum = cumsum(CDD20))

sites_climate_list.df.cumsum$CDD20.cumsum %>% range


########################
#Define growing seasons#
sites_climate_list.df.cumsum$season = vector(mode = "character", length = length(sites_climate_list.df.cumsum$tmin))

sites_climate_list.df.cumsum[,colnames(sites_climate_list.df.cumsum) == "season"][
    sites_climate_list.df.cumsum[,colnames(sites_climate_list.df.cumsum) == "HDD4.cumsum"] < 100
] = "nongrowing"

sites_climate_list.df.cumsum[,colnames(sites_climate_list.df.cumsum) == "season"][
    sites_climate_list.df.cumsum[,colnames(sites_climate_list.df.cumsum) == "HDD4.cumsum"] >= 100 &
    sites_climate_list.df.cumsum[,colnames(sites_climate_list.df.cumsum) == "CDD20.cumsum"] > -500
] = "growing"

sites_climate_list.df.cumsum[,colnames(sites_climate_list.df.cumsum) == "season"][
    sites_climate_list.df.cumsum[,colnames(sites_climate_list.df.cumsum) == "CDD20.cumsum"] <= -500
] = "nongrowing"

################
#Calc HDD4, freezeThaw and ppt for seasons

sites_climate.seasonal_summaries = sites_climate_list.df.cumsum %>%
    group_by(year, Site, season) %>%
    summarize(HDD4.sum = sum(HDD4), freezeThaw.sum = sum(freezeThaw), ppt.sum = sum(ppt))

#get min of HDD4.cumsum and min (bc negative) of CDD20.cumsum
#will then do filtering join based on extracted vals
temps_for_season_start = sites_climate_list.df.cumsum %>%
    filter(season == "growing") %>%
    group_by(year, Site) %>%
    summarize(HDD4.cumsum = min(HDD4.cumsum))
temps_for_season_end = sites_climate_list.df.cumsum %>%
    filter(season == "growing") %>%
    group_by(year, Site) %>%
    summarize(CDD20.cumsum = min(CDD20.cumsum))
year_lens = sites_climate_list.df.cumsum %>%
    group_by(year, Site) %>%
    summarize(year_len = max(doy) %>% as.numeric)

#Growing season start and end dates and length
grow_start = semi_join(sites_climate_list.df.cumsum, temps_for_season_start) %>%
    group_by(Site, year) %>%
    summarize(doy.grow_start = min(doy) %>% as.numeric)
grow_end = semi_join(sites_climate_list.df.cumsum, temps_for_season_end) %>%
    group_by(Site, year) %>%
    summarize(doy.grow_end = min(doy) %>% as.numeric)

grow_season_dates = full_join(grow_start, grow_end) %>%
    full_join(., year_lens)
grow_season_dates$length.growing = grow_season_dates$doy.grow_end - grow_season_dates$doy.grow_start + 1
grow_season_dates$length.nongrowing = grow_season_dates$year_len - grow_season_dates$length.growing

#pivot_longer season
grow_season_dates.long = grow_season_dates %>%
    pivot_longer(cols = c(length.growing, length.nongrowing), names_to = "season", names_prefix = "length.", values_to = "length")

#join climate with season lens
sites_climate.seasonal_summaries = full_join(sites_climate.seasonal_summaries, grow_season_dates.long)

#order factor
#sites_climate.seasonal_summaries$Site = factor(sites_climate.seasonal_summaries$Site, levels = c("MEN1", "MES1", "ADN1", "ADS1", "CW1", "TSP1", "GK1j", "ASH2", "MI1", "WF1"))

#means and sd
sites_climate.seasonal_summaries.means = sites_climate.seasonal_summaries %>%
    group_by(Site, season) %>%
    summarize(
        HDD4.mean = mean(HDD4.sum),
        freezeThaw.mean = mean(freezeThaw.sum),
        HDD4.sd = sd(HDD4.sum),
        freezeThaw.sd = sd(freezeThaw.sum),
        ppt.mean = mean(ppt.sum),
        ppt.sd = sd(ppt.sum),
        length.mean = mean(length),
        length.sd = sd(length)

)

saveRDS(sites_climate.seasonal_summaries, file = "sample_metadata/summarized_12_year_climate_dailys.rds")
saveRDS(sites_climate.seasonal_summaries.means, file = "sample_metadata/summarized_12_year_climate_dailys.means.rds")

write.table(sites_climate.seasonal_summaries, file = "sample_metadata/summarized_12_year_climate_dailys.txt", quote = F, row.names = F, sep = "\t")
write.table(sites_climate.seasonal_summaries.means, file = "sample_metadata/summarized_12_year_climate_dailys.means.txt", quote = F, row.names = F, sep = "\t")

sites_climate.seasonal_summaries.means.wide = sites_climate.seasonal_summaries.means %>% pivot_wider(names_from = season, values_from = c(-Site, -season))

write.table(sites_climate.seasonal_summaries.means.wide, file = "sample_metadata/site_cliamte.GDD.txt", quote = F, row.names = F, sep = "\t")

#################
#Plots

##########
#GDD4

p1 = ggplot(sites_climate.seasonal_summaries, aes(year, HDD4.sum, color = Site)) +
geom_line() +
geom_point() +
facet_wrap(~season, nrow = 2, scales = "free") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")) +
my_gg_theme +
labs(
    y = expression(paste("Cumulative GDD"[4])),
    title = "Growing season vs. non-growing season\ncumulative GDD"
)

pdf("PRISM_maps/GDD4_over_time.pdf", width = 10, height = 6)
p1
dev.off()

p1 = ggplot(sites_climate.seasonal_summaries.means, aes(Site, HDD4.mean, color = Site)) +
geom_errorbar(aes(ymin = HDD4.mean-HDD4.sd, ymax = HDD4.mean+HDD4.sd), color = "black", width = 0.3) +
geom_point(size = 3) +
facet_wrap(~season, nrow = 2, scales = "free") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")) +
my_gg_theme +
scale_x_discrete(labels = NULL) +
labs(
y = expression(paste("Cumulative GDD"[4])),
title = "Growing season vs. non-growing season\ncumulative GDD (12 yr mean +/- sd)"
)

pdf("PRISM_maps/GDD4_12_yr_mean.pdf", width = 10, height = 6)
p1
dev.off()

#############
#Freeze-thaw

p1 = ggplot(sites_climate.seasonal_summaries, aes(year, freezeThaw.sum, color = Site)) +
geom_line() +
geom_point() +
facet_wrap(~season, nrow = 2, scales = "free") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")) +
my_gg_theme +
labs(
y = expression(paste("Single day freeze-thaw")),
title = "Growing season vs. non-growing season\nsingle day freeze-thaw events"
)

pdf("PRISM_maps/freeze-thaw_over_time.pdf", width = 10, height = 6)
p1
dev.off()

p1 = ggplot(sites_climate.seasonal_summaries.means, aes(Site, freezeThaw.mean, color = Site)) +
geom_errorbar(aes(ymin = freezeThaw.mean-freezeThaw.sd, ymax = freezeThaw.mean+freezeThaw.sd), color = "black", width = 0.3) +
geom_point(size = 3) +
facet_wrap(~season, nrow = 2, scales = "free") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")) +
my_gg_theme +
scale_x_discrete(labels = NULL) +
labs(
y = expression(paste("Cumulative GDD"[4])),
title = "Growing season vs. non-growing season\nsingle day freeze-thaw events (12 yr mean +/- sd)"
)

pdf("PRISM_maps/freeze-thaw_12_yr_mean.pdf", width = 10, height = 6)
p1
dev.off()

######
#ppt

p1 = ggplot(sites_climate.seasonal_summaries, aes(year, ppt.sum, color = Site)) +
geom_line() +
geom_point() +
facet_wrap(~season, nrow = 2, scales = "free") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")) +
my_gg_theme +
labs(
y = "Precipitation (mm)",
title = "Growing season vs. non-growing season\nprecipitation"
)

pdf("PRISM_maps/ppt_over_time.pdf", width = 10, height = 6)
p1
dev.off()

p1 = ggplot(sites_climate.seasonal_summaries.means, aes(Site, ppt.mean, color = Site, shape = season)) +
geom_errorbar(aes(ymin = ppt.mean-ppt.sd, ymax = ppt.mean+ppt.sd), color = "black", width = 0.3, position = position_dodge(width = 0.9)) +
geom_point(size = 3, position = position_dodge(width = 0.9)) +
#p1 = ggplot(sites_climate.seasonal_summaries.means, aes(Site, ppt.mean, color = Site, fill = season)) +
#geom_boxplot(position = position_dodge(width = 0.9)) +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")) +
my_gg_theme +
scale_fill_manual(values = c("growing" = "white", "nongrowing" = "light grey")) +
scale_x_discrete(labels = NULL) +
labs(
y = "Precipitation (mm)",
title = "Growing season vs. non-growing season\ntotal precipitation (12 yr mean +/- sd)"
)

pdf("PRISM_maps/ppt_12_yr_avg.pdf", width = 10, height = 6)
p1
dev.off()


##################
#Growing season len

p1 = ggplot(sites_climate.seasonal_summaries %>% filter(season == "growing"), aes(year, length, color = Site)) +
geom_line() +
geom_point() +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")) +
my_gg_theme +
labs(
y = expression(paste("Growing season length (days)")),
title = "Growing season length"
)


pdf("PRISM_maps/growing_season_length_over_time.pdf", width = 10, height = 5)
p1
dev.off()


p1 = ggplot(sites_climate.seasonal_summaries.means %>% filter(season == "growing"), aes(Site, length.mean, color = Site)) +
geom_errorbar(aes(ymin = length.mean-length.sd, ymax = length.mean+length.sd), color = "black", width = 0.3) +
geom_point(size = 3) +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")) +
my_gg_theme +
scale_x_discrete(labels = NULL) +
labs(
y = "Growing season length (days)",
title = "Growing season length (12 yr mean +/- sd)"
)

pdf("PRISM_maps/growing_season_length_12_yr_mean.pdf", width = 10, height = 5)
p1
dev.off()

sites_climate.season_sum_w_len


####################
#Grow season len normalized GDD

p1 = ggplot(sites_climate.seasonal_summaries, aes(year, HDD4.sum/length, color = Site)) +
geom_point() +
geom_line() +
facet_wrap(~season, ncol = 1, scales = "free") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")) +
my_gg_theme +
scale_x_discrete(labels = NULL) +
labs(
y = expression(paste("GDD"[4], " day"^-1)),
title = "Daily mean GDD"
)

pdf("PRISM_maps/growing_season_daily_mean_GDD_over_time.pdf", width = 10, height = 6)
p1
dev.off()

p1 = ggplot(sites_climate.seasonal_summaries, aes(Site, HDD4.sum/length, color = Site)) +
geom_boxplot() +
facet_wrap(~season, ncol = 1, scales = "free") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")) +
my_gg_theme +
scale_x_discrete(labels = NULL) +
labs(
y = expression(paste("GDD"[4], " day"^-1)),
title = "Daily GDD 12 yr average"
)

pdf("PRISM_maps/growing_season_daily_GDD_12_yr.pdf", width = 10, height = 6)
p1
dev.off()


####################
#Grow season len normalized ppt

p1 = ggplot(sites_climate.seasonal_summaries, aes(year, ppt.sum/length, color = Site)) +
geom_point() +
geom_line() +
facet_wrap(~season, ncol = 1, scales = "free") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")) +
my_gg_theme +
scale_x_discrete(labels = NULL) +
labs(
y = expression(paste("precipitation mm day"^-1)),
title = "Daily precipitation"
)

pdf("PRISM_maps/growing_season_daily_ppt_over_time.pdf", width = 10, height = 6)
p1
dev.off()

p1 = ggplot(sites_climate.seasonal_summaries.means, aes(Site, ppt.mean/length.mean, color = Site, shape = season)) +
#geom_boxplot(position = position_dodge(width = 0.9)) +
geom_errorbar(aes(
    #the term for combining errors (square root of sum of squared proportional error)
        ymin = ppt.mean/length.mean - sqrt((ppt.sd/ppt.mean)^2+(length.sd/length.mean)^2),
        ymax = ppt.mean/length.mean + sqrt((ppt.sd/ppt.mean)^2+(length.sd/length.mean)^2)
    ),
    color = "black",
    position = position_dodge(width = 0.9),
    width = 0.3
) +
geom_point(position = position_dodge(width = 0.9), size = 3) +
#facet_wrap(~season, ncol = 1, scales = "free") +
scale_color_brewer(palette = "Paired", name = "Site", labels = c("ADS1" = "NY south", "ADN1" ="NY north", "MEN1" = "ME north", "MES1" = "ME south", "CW1" = "NH", "TSP1" = "PA", "ASH2" = "NC", "GK1j" = "WV", "MI1" = "MI", "WF1" = "WI")) +
my_gg_theme +
scale_x_discrete(labels = NULL) +
scale_fill_manual(values = c("growing" = "white", "nongrowing" = "light grey")) +
labs(
y = expression(paste("precipitation mm day"^-1)),
title = "Growing season vs. non-growing season\ndaily precipiation 12 yr mean +/- sd"
)

pdf("PRISM_maps/growing_season_mean_daily_ppt_12_yr.pdf", width = 10, height = 6)
p1
dev.off()
