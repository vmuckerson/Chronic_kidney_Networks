library(cluster)
library(dbscan)
library(BiocManager)
library(cluster)
library(mclust)
library(tidyverse)
library(ConsensusClusterPlus)
library(ALL)
library(dplyr)

all_data <- readRDS(file = "Documents/all_data.RDS")
only_num <- readRDS(file = "Documents/only_num.RDS")
comdata2 <- readRDS(file = "Documents/comdata2.RDS")
all_pca <- readRDS(file = "Documents/all_pca.RDS")
com_pca2 <- readRDS(file = "Documents/com_pca2.RDS")
dataScaled <- scale(only_num[,])
disease <- all_data[,4]
experiment <- all_data[,1]
platform <- all_data[,2]
comScaled <- scale(t(comdata2))

#K-MEANS CLUSTERING ---
#CHOOSING K -----------
k <- list()
for(i in 1:15){
  k[[i]] <-kmeans(na.omit(dataScaled[,-1]), i)
}

k

betweenss_totss <- list()
for(i in 1:15){
  betweenss_totss[[i]] <- k[[i]]$betweenss/k[[i]]$totss
}

plot(1:15, betweenss_totss, type = "b",
     ylab = "Between SS / Total SS", xlab = "Number of Clusters (k)")




#CHOOSING K POST COMBAT -----------
k2 <- list()
for(i in 1:15){
  k2[[i]] <-kmeans(comScaled[,], i)
}

k2

betweenss_totss2 <- list()
for(i in 1:15){
  betweenss_totss2[[i]] <- k2[[i]]$betweenss/k2[[i]]$totss
}

plot(1:15, betweenss_totss2, type = "b",
     ylab = "Between SS / Total SS", xlab = "Number of Clusters (k)")

#CLUSTERING -------------
fitK <- kmeans(dataScaled[,-1], 5)

tmpk <- data.frame(PC1= all_pca$x[,1], PC2=all_pca$x[,2],clust=fitK$cluster, dis=disease, ex = experiment, plat = platform)
ggplot(tmpk, aes(PC1,PC2,col=as.factor(disease))) +
  geom_point(aes(shape=as.factor(fitK$cluster)))+
  stat_ellipse(aes(PC1, PC2, color = as.factor(fitK$cluster)), type = "norm")


## CLUSTERING POST COMBAT--------
fitKc <- kmeans(comScaled[,], 3)

tmpkc <- data.frame(PC1= com_pca2$x[,1], PC2=com_pca2$x[,2],clust=fitKc$cluster, dis=disease, ex = experiment, plat = platform)
ggplot(tmpkc, aes(PC1, PC2, col=as.factor(disease)))+
  geom_point(aes(shape=as.factor(fitKc$cluster))) +
  stat_ellipse(aes(PC1, PC2, color = as.factor(fitKc$cluster)), type = "norm")


#HIERARCHICAL CLUSTERING --------
#PRE COMBAT DATA ----------------
d <- dist(dataScaled[,-1])
fitH <- hclust(d)
plot(fitH)
rect.hclust(fitH, k = 5, border = "red")
clusters <- cutree(fitH, k = 5)


tmp <- data.frame(PC1= all_pca$x[,1], PC2=all_pca$x[,2],clust=clusters, dis=disease, ex = experiment, plat = platform)
ggplot(tmp, aes(PC1,PC2,col=as.factor(disease))) +
  geom_point(aes(shape=as.factor(clust))) +
  stat_ellipse(aes(PC1, PC2, color = as.factor(clust)), type = "norm")



#fitH2 <- hclust(d, "ward.D2")
#plot(fitH2)
#rect.hclust(fitH2, k = 4, border = "red")
#clusters2 <- cutree(fitH2, k = 4)

#tmp2 <- data.frame(pc.x= all_pca$x[,1], pc.y=all_pca$x[,2],clust=clusters2, dis=disease)
#ggplot(tmp2, aes(pc.x,pc.y,col=as.factor(disease))) + geom_point(aes(shape=as.factor(clust)))



#POST COMBAT DATA--------------
dc <- dist(comScaled)
fitHc <- hclust(dc)
plot(fitHc)
rect.hclust(fitHc, k = 3, border = "red")
clustersc <- cutree(fitHc, k = 3)
plot(comScaled, col = clustersc)


tmpc <- data.frame(PC1= com_pca2$x[,1], PC2= com_pca2$x[,2],clust=clustersc, dis=disease, ex = experiment, plat = platform)
ggplot(tmpc, aes(PC1,PC2,col=as.factor(disease))) +
  geom_point(aes(shape=as.factor(clust))) +
  stat_ellipse(aes(PC1, PC2, color = as.factor(clust), type = "norm"))



#fitH2c <- hclust(dc, "ward.D2")
#plot(fitH2c)
#rect.hclust(fitH2c, k = 4, border = "red")
#clusters2c <- cutree(fitH2c, k = 4)
#plot(comdata, col = clusters2c)

#tmpc2 <- data.frame(pc.x= com_pca2$x[,1], pc.y=com_pca2$x[,2],clust=clusters2c, dis=disease, ex = experiment, plat = platform)
#ggplot(tmpc2, aes(pc.x,pc.y,col=as.factor(experiment))) +
#  geom_point(aes(shape=as.factor(clust)))


#DENSITY BASED CLUSTERING -------
kNNdistplot(dataScaled, k = 5)
abline(h = 52, col = "red", lty = 2)
fitD <- dbscan(dataScaled, eps = 52, minPts = 5)
fitD
tmpd <- data.frame(PC1= all_pca$x[,1], PC2=all_pca$x[,2],clust=fitD$cluster, dis=disease, ex = experiment, plat = platform)
ggplot(tmpd, aes(PC1,PC2,col=as.factor(disease))) +
  geom_point(aes(shape=as.factor(fitD$cluster)))


#DBC ON POSTCOMBAT DATA---------
kNNdistplot(comScaled, k = 3)
abline(h = 85, col = "red", lty = 2)
fitDc <- dbscan(comScaled, eps = 85, minPts = 3)
fitDc
tmpdc <- data.frame(PC1= com_pca2$x[,1], PC2=com_pca2$x[,2],clust=fitDc$cluster, dis=disease, ex = experiment, plat = platform)
ggplot(tmpdc, aes(PC1,PC2,col=as.factor(disease))) +
  geom_point(aes(shape=as.factor(fitDc$cluster)))



#MODEL-BASED CLUSTERING ---------
fitM <- Mclust(dataScaled[,-1])
plot(fitM)
tmpM <- data.frame(PC1= com_pca2$x[,1], PC2=com_pca2$x[,2],clust=fitM$classification, dis=disease, ex = experiment, plat = platform)
ggplot(tmpM, aes(PC1,PC2,col=as.factor(disease))) +
  geom_point(aes(shape=as.factor(fitM)))

#MBC ON POST COMBAT DATA---------
fitM2 <- Mclust(comScaled)
plot(fitM2)
tmpMc <- data.frame(PC1= com_pca2$x[,1], PC2=com_pca2$x[,2],clust=fitM2$classification, dis=disease, ex = experiment, plat = platform)
ggplot(tmpMc, aes(PC1,PC2,col=as.factor(disease))) +
  geom_point(aes(shape=as.factor(fitM2$classification)))

#Consensus Clustering------------
clu <- t(dataScaled)
mads = apply(clu,1,mad)
d=clu[rev(order(mads))[1:5000],]
d=sweep(d,1,apply(d,1,median,na.rm=T))
title=tempdir()
results=ConsensusClusterPlus(d, maxK=12,
                             reps=1000,
                             pItem = 0.8,
                             pFeature=1,
                             title="Consensus Cluster on Precombat Data",
                             innerLinkage="average",
                             finalLinkage = "average",
                             clusterAlg="hc",
                             distance="pearson")


icl <- calcICL(results)


#CONSENSUS CLUSTERING ON POST COMBAT DATA------
cluc <- t(comScaled)
mads = apply(cluc, 1, mad)
dc=cluc[rev(order(mads))[1:5000],]
dc=sweep(dc,1,apply(dc,1,median,na.rm=T))
resultsc=ConsensusClusterPlus(dc, maxK=12,
                             reps=1000,
                             pItem = 0.8,
                             pFeature=1,
                             title="Consensus Cluster on Postcombat Data",
                             innerLinkage="average",
                             finalLinkage = "average",
                             clusterAlg="hc",
                             distance="pearson")

iclc = calcICL(resultsc)

