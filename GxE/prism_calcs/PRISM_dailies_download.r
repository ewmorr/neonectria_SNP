require(prism)

prism_set_dl_dir("/mnt/home/garnas/ericm/PRISM_dailies_tmin")
get_prism_dailys(type = "tmin", minDate="2007-01-01", maxDate="2018-12-31", keepZip = TRUE)
prism_set_dl_dir("/mnt/home/garnas/ericm/PRISM_dailies_tmax")
get_prism_dailys(type = "tmax", minDate="2007-01-01", maxDate="2018-12-31", keepZip = TRUE)

