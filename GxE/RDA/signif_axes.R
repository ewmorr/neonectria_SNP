library(vegan)

args = commandArgs(trailingOnly=TRUE)
Y.rda = readRDS(args[1])

start.time <- Sys.time()
signif.axis <- anova.cca(Y.rda, by="axis", parallel=getOption("mc.cores"))
end.time <- Sys.time()

saveRDS( file.path(dirname(args[1]), "signif_axis.rds"))
