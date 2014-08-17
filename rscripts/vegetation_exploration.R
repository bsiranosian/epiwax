setwd("~/GitHub/epiwax")
library(robCompositions)

# load vegetation component distribution data
veg.fa <- as.matrix(read.table('distribution_data//veg_acid.txt',header = T,check.names = F))
veg.fal <- as.matrix(read.table('distribution_data//veg_alcohol.txt',header = T,check.names = F))
veg.alk <- as.matrix(read.table('distribution_data//veg_alkane.txt',header = T,check.names = F))

#subset each to the interesting components because we don't have enough samples
veg.fa <- veg.fa[,c(1,3,5,7,9,11,13)]

#normalize each by rowsums
veg.fa <- veg.fa/rowSums(veg.fa)
veg.fal <- veg.fal/rowSums(veg.fal)
veg.alk <- veg.alk/rowSums(veg.alk)

veg.fa.p1 <- pcaCoDa(veg.fa)
