# trying robust PCA on compositional data
# percent compoistion of varying component classes in wax samples. 

setwd("~/GitHub/epiwax")
library(robCompositions)
# load FAL
fal.conc <- as.matrix(read.table('distribution_data//transect_alcohol.txt',header = T,check.names = F))
fal.even <- fal.conc/rowSums(fal.conc)

#load FA
fa.conc <- as.matrix(read.table('distribution_data//transect_acid.txt', header=T,check.names = F))
fa.conc.even <- fa.conc[,c(1,3,5,7,9,11,13,15)]
fa.even <- fa.conc.even/rowSums(fa.conc.even)

#load ALK
alk.conc <- as.matrix(read.table('distribution_data//transect_alkane.txt', header=T,check.names = F))
alk.conc.odd <- alk.conc[,c(1,3,5,7,9)]
alk.odd <- alk.conc.odd/rowSums(alk.conc.odd)

#plot some distributions
# FAL EVEN
pdf('figures//component_distributions/transect_alcohol.pdf',width=8.5, height=11)
par(mfrow=c(5,4), mai = c(0.45,0.45,0.45,0.45), mar=c(2,2,2,2))
for (i in seq(1,20)){
  barplot(fal.even[i,], main = i)
}
dev.off()
# FA EVEN
pdf('figures//component_distributions/transect_acid.pdf',width=8.5, height=11)
par(mfrow=c(5,4), mai = c(0.45,0.45,0.45,0.45), mar=c(2,2,2,2))
for (i in seq(1,20)){
  barplot(fa.even[i,], main = i)
}
dev.off()
# ALK ODD
pdf('figures//component_distributions/transect_alkane.pdf',width=8.5, height=11)
par(mfrow=c(5,4), mai = c(0.45,0.45,0.45,0.45), mar=c(2,2,2,2))
for (i in seq(1,20)){
  barplot(alk.odd[i,], main = i)
}
dev.off()

# Distance between samples, heatmaps
par(mfrow=c(2,2))
heatmap(as.matrix(dist(fal.even)),main = "alcohol")
heatmap(as.matrix(dist(fa.even)), main= "acid")
heatmap(as.matrix(dist(alk.odd)), main= "alkane")


# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name
d <- dist(fa.even) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=8) # k is the number of dim
fit # view results
# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric  MDS",	type="n")
text(x, y, labels = row.names(fa.even), cex=.7)

#### starting to try out robCompositions 
ternaryDiag(fal[,1:3])

#pca alcohol
fal.even.p1 <- pcaCoDa(fal.even)
plot(fal.even.p1)

# fatty acid
fa.even.p1 <- pcaCoDa(fa.even)
plot(fa.even.p1)

alk.odd.p1 <- pcaCoDa(alk.odd)
plot(alk.odd.p1)
alk.odd.p1
