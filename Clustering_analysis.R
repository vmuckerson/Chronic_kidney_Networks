if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("cluster")
library(BiocManager)
library(cluster)


dataScaled <- scale(only_num[,])

#K-MEANS CLUSTERING ---
#CHOOSING K -----------
k <- list()
for(i in 1:15){
  k[[i]] <-kmeans(only_num[,], i)
}

k

betweenss_totss <- list()
for(i in 1:15){
  betweenss_totss[[i]] <- k[[i]]$betweenss/k[[i]]$totss
}

plot(1:15, betweenss_totss, type = "b",
     ylab = "Between SS / Total SS", xlab = "Number of Clusters (k)")


#CHOOSING K POST COMBAT -----------
comScaled <- scale(comdata)

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
fitK <- kmeans(dataScaled[,], 2)
plot(only_num, col = fitK$cluster)

fitK <- kmeans(dataScaled[,], 4)
plot(only_num, col = fitK$cluster)

fitK <- kmeans(dataScaled[,], 6)
plot(only_num, col = fitK$cluster)

fitK <- kmeans(all.pca[,], 7)
plot(only_num, col = fitK$cluster)

## CLUSTERING POST COMBAT--------
fitK2 <- kmeans(comScaled[,], 2)
plot(comdata, col = fitK2$cluster)

fitK2 <- kmeans(comScaled[,], 4)
plot(comdata, col = fitK2$cluster)

fitK2 <- kmeans(comScaled[,], 6)
plot(comdata, col = fitK2$cluster)

fitK2 <- kmeans(comScaled[,], 7)
plot(comdata, col = fitK2$cluster)



#HIERARCHICAL CLUSTERING --------






#MODEL-BASED CLUSTERING ---------
library(mclust)
fitM <- Mclust(only_num)
plot(fitM)

#MBC ON POST COMBAT DATA---------
fitM2 <- Mclust(comdata)
plot(fitM2)









