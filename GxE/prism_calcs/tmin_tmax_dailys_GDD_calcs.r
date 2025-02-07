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

write.table(sites_climate.seasonal_summaries.means.wide, file = "sample_metadata/site_climate.GDD.txt", quote = F, row.names = F, sep = "\t")
