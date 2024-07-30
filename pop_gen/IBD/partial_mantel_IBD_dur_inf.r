require(vegan)
require(dplyr)
require(reshape2)
require(ggplot2)
source("R_scripts/ggplot_theme.txt")


site_info = read.csv("data/sample_metadata/site_info.dur_inf.csv")

#subset and order
sites = c("CT",
"ME.S", "NC", "NH.BART", "NH.CW", "NY.N", "NY.S", "NY.W", "PA",
"PA.W", "VA", "WV")
dur_inf.subset = site_info %>% filter(state.name %in% sites) %>% select(state.name, duration_infection)
#range and mean of durInf without VA
dur_inf.subset$duration_infection[-9] %>% range
dur_inf.subset$duration_infection[-9] %>% mean


#need to order the df to compare distances correctly
dur_inf.subset.order = dur_inf.subset[match(sites, dur_inf.subset$state.name), ]
rownames(dur_inf.subset.order) = dur_inf.subset.order$state.name

dur_inf.dist = dist(dur_inf.subset.order$duration_infection, method = "euclidean", diag = FALSE, upper = FALSE)

#import gen and geo distance
distances.list = readRDS("data/intermediate_RDS/DSCE.five_samples_per_site.rds")
mean_Dgen = Reduce("+", distances.list) / length(distances.list)

Dgeo = readRDS("data/intermediate_RDS/geo_distance.min_5_samples.rds")


#mean_Dgen = structure(c(0.327054284813512, 0.33019388260635, 0.311480420341884,
0.320345373087978, 0.327170825312449, 0.318511606142414, 0.310771628900295,
0.329864438413669, 0.321813865417881, 0.374568119282915, 0.340920031009168,
0.32740606392453, 0.314750048967457, 0.321121989441771, 0.323546104641977,
0.316927942315207, 0.326976764055849, 0.330524373664444, 0.32812186851919,
0.378302148757357, 0.339948048565828, 0.322491891962815, 0.328781231479944,
0.332708692782155, 0.325569362424847, 0.329881224246572, 0.33936749015364,
0.332237584913339, 0.374875187719548, 0.335695730863176, 0.313830585735777,
0.317192693142494, 0.308575339871422, 0.307711006603123, 0.323979680607829,
0.313837184108629, 0.369042013095151, 0.331632380708735, 0.321969915435215,
0.313652683275849, 0.319878668418307, 0.328463493990181, 0.321872578681901,
0.374929330379957, 0.33779794745129, 0.310953612841122, 0.326229911701076,
0.326564068536057, 0.322982280443371, 0.38078579047822, 0.339339588324459,
0.316188144117529, 0.31948739776192, 0.31582903345406, 0.37321359705047,
0.331711843141456, 0.328034865396458, 0.315666345800262, 0.373610851468843,
0.339688307960393, 0.329248151117817, 0.383354386138072, 0.346124276373405,
0.375008971370167, 0.339392151826257, 0.373290216667519), Labels = c("CT",
"ME.S", "NC", "NH.BART", "NH.CW", "NY.N", "NY.S", "NY.W", "PA",
"PA.W", "VA", "WV"), Size = 12L, class = "dist", Diag = FALSE, Upper = FALSE, method = "Edwards")

#Dgeo = structure(c(647867.64331137, 1459093.77642383, 359233.889326968,
278944.401616104, 393611.221841495, 227982.167564933, 398327.059881948,
295694.785989413, 657332.64055266, 816245.79577908, 915018.076746268,
2106955.61116748, 323712.25040289, 370719.975123838, 606761.406619077,
704096.511567019, 953236.578051644, 934751.083545949, 1225840.71616558,
1464087.78917782, 1559011.88689922, 1804528.81619937, 1737276.06000047,
1656384.35769276, 1457191.63598701, 1240678.08494351, 1178772.8489121,
1043511.2093198, 642891.180646641, 566880.653369328, 146811.149326077,
311632.863183784, 380421.302011604, 629564.548523671, 625961.246379126,
902168.256876964, 1162992.42660029, 1248543.47514078, 401693.3328347,
388536.470766616, 626721.293027877, 571500.629812364, 899648.158889279,
1094623.2074738, 1193613.13585933, 221949.243778983, 417735.024909979,
520252.979315865, 664781.90146369, 1036503.10582002, 1089511.15453256,
249143.919731233, 299737.786616343, 522243.290019075, 826753.989289335,
892401.525736529, 221910.884051497, 274078.226591962, 634743.178350607,
674111.30935297, 402976.819030762, 538688.67312651, 624342.881849541,
519738.887021275, 498924.428982888, 150205.297652232), Size = 12L, Labels = c("CT",
"ME.S", "NC", "NH.BART", "NH.CW", "NY.N", "NY.S", "NY.W", "PA",
"PA.W", "VA", "WV"), Diag = FALSE, Upper = FALSE, method = "euclidean", class = "dist")

mantel.partial(mean_Dgen, Dgeo, dur_inf.dist)
mantel.partial(mean_Dgen, log(Dgeo), dur_inf.dist)

mantel(mean_Dgen, Dgeo)
mantel(mean_Dgen, log(Dgeo))

mantel(mean_Dgen, dur_inf.dist)
#r = 0.4647, P = 0.002

mantel(Dgeo, dur_inf.dist)
#r = 0.7147, P = 0.001

mantel.partial(mean_Dgen, dur_inf.dist, Dgeo)
mantel.partial(mean_Dgen, dur_inf.dist, log(Dgeo))



Dgeo.long = subset(reshape2::melt(Dgeo %>% as.matrix), value != 0)
Dgen.long = subset(reshape2::melt(mean_Dgen %>% as.matrix), value != 0)
Dinf.long = subset(reshape2::melt(dur_inf.dist %>% as.matrix), value != 0)

Dgen.Dgeo = data.frame(gen = Dgen.long$value, geo = Dgeo.long$value, inf = Dinf.long$value)

p1 = ggplot(Dgen.Dgeo, aes(x = inf, y = gen)) +
geom_point() +
geom_smooth(method = "lm", color = "black", se = F, linetype = 2) +
labs(x = "Pairwise site duration infection distance (years)", y = expression(paste("Pairwise site genetic distance (D"["CSE"], ")"))) +
my_gg_theme

p2 = ggplot(Dgen.Dgeo, aes(x = geo/1000, y = inf)) +
geom_point() +
geom_smooth(method = "lm", color = "black", se = F, linetype = 2) +
labs(x = "Pairwise site geographic distance (km)", y = "Pairwise site duration\ninfection distance (years)") +
#scale_x_continuous(trans = "log", breaks = c(200, 500, 1000, 2000)) +
my_gg_theme


pdf("figures/IBT.Dcse_v_dur_inf.pdf", width = 9, height = 5)
p1
dev.off()

pdf("figures/IBT.Dgeo_v_dur_inf.pdf", width = 9, height = 5)
p2
dev.off()
