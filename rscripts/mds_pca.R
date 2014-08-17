# trying robust PCA on compositional data
# percent compoistion of varying component classes in wax samples. 

setwd("~/GitHub/epiwax")
library(robCompositions)
# load FAL
fal_conc <- as.matrix(read.table('distribution_data//fatty_alcohol.txt',header = T))
falt <- apply(fal_conc, 2, function(x) x/sum(x))
fal <- t(falt)

#load FA
fa_conc <- as.matrix(read.table('distribution_data//fatty_acid.txt', header=T))
fat <- apply(fa_conc, 2, function(x) x/sum(x))
fa <- t(fat)
fa_even <- fa[,colnames(fa)]

#load ALK
alk_conc <- as.matrix(read.table('distribution_data//alkane.txt', header=T))
alkt <- apply(alk_conc, 2, function(x) x/sum(x))
alk <- t(alkt)


# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name
d <- dist(fal) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results
# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric  MDS",	type="n")
text(x, y, labels = row.names(t(fal)), cex=.7)

#### starting to try out robCompositions 
ternaryDiag(fal[,1:3])

#pca alcohol
fal_p1 <- pcaCoDa(fal)
plot(fal_p1)

# fatty acid
fa_p1 <- pcaCoDa(fa)
plot(fa_p1)
# too small sample size
fa.subset <- fa[,as.integer(colnames(fa)) %in% c(20,22,24,26,28,30,32,34)]
fa.subset_p1 <- pcaCoDa(fa.subset)
plot(fa.subset_p1)

alk_p1 <- pcaCoDa(alk)
# too small sample size? also some zeros on evens
alk.subset <- alk[, as.integer(colnames(alk)) %in% c(25,27,29,31)]
alk.subset_p1 <- pcaCoDa(alk.subset)
plot(alk.subset_p1)
