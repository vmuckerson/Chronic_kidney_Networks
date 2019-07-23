library(cluster)
library(dbscan)
library(BiocManager)
library(cluster)
library(mclust)
library(tidyverse)
library(ConsensusClusterPlus)
library(ALL)
library(dplyr)
library(ComplexHeatmap)
library(formattable)
library(ggbiplot)

setwd("/Users/saezlab/Documents/CKD_Data/Glom_Data/")

#Load data
data1 <- readRDS(file = "./all_104948.rds")
pca1 <- readRDS(file = "./pca_104948.rds")
Pdata1 <- readRDS(file = "./Pdata_104948.rds")
scdata1 <- scale(data1)

data2 <- readRDS(file = "./all_GSE37460.rds")
pca2 <- readRDS(file = "./pca_GSE37460.rds")
Pdata2 <- readRDS(file = "./Pdata_GSE37460.rds")
scdata2 <- scale(data2)

data3 <- readRDS(file = "./all_GSE32591.rds")
pca3 <- readRDS(file = "./pca_GSE32591.rds")
Pdata3 <- readRDS(file = "./Pdata_GSE32591.rds")
scdata3 <- scale(data3)


#Selecting most variable genes
var1 <- rev(sort(apply(data1, 1, FUN = var)))
var2 <- rev(sort(apply(data2, 1, FUN = var)))
var3 <- rev(sort(apply(data3, 1, FUN = var)))

top1 <- var1[1:2000]
top2 <- var2[1:2000]
top3 <- var3[1:2000]

data1 <- scale(data1[rownames(data1) %in% names(top1),])
data2 <- scale(data2[rownames(data2) %in% names(top2),])
data3 <- scale(data3[rownames(data3) %in% names(top3),])



##### CLUSTERING #############

#Consensus Clustering------------

mad1 <- apply(data1, 1, mad)
d1 <- data1[rev(order(mad1))[1:1000],]
d1 <- sweep(d1,1,apply(d1,1,median,na.rm=TRUE))
d1 <- na.omit(d1)
results1 <- ConsensusClusterPlus(d1, maxK=6,
                                 reps=1000,
                                 pItem = 0.8,
                                 pFeature=1,
                                 title="Consensus Cluster on GSE104948 most variable genes",
                                 innerLinkage="average",
                                 finalLinkage = "average",
                                 clusterAlg="hc",
                                 distance="pearson")
icl1 <- calcICL(results1)



mad2 <- apply(data2, 1, mad)
d2 <- data2[rev(order(mad2))[1:5000],]
d2 <- sweep(d2,1,apply(d2,1,median,na.rm=TRUE))
d2 <- na.omit(d2)
results2 <- ConsensusClusterPlus(d2, maxK=6,
                                 reps=1000,
                                 pItem = 0.8,
                                 pFeature=1,
                                 title="Consensus Cluster on GSE37460 most variable genes",
                                 innerLinkage="average",
                                 finalLinkage = "average",
                                 clusterAlg="hc",
                                 distance="pearson")
icl2 <- calcICL(results2)



mad3 <- apply(data3, 1, mad)
d3 <- data3[rev(order(mad3))[1:5000],]
d3 <- sweep(d3,1,apply(d3,1,median,na.rm=TRUE))
d3 <- na.omit(d3)
results3 <- ConsensusClusterPlus(d3, maxK=6,
                                 reps=1000,
                                 pItem = 0.8,
                                 pFeature=1,
                                 title="Consensus Cluster on GSE32591 most variable genes",
                                 innerLinkage="average",
                                 finalLinkage = "average",
                                 clusterAlg="hc",
                                 distance="pearson")
icl3 <- calcICL(results3)


##### HEAT MAPS #########

#CLUSTERING HEAT MAPS------------------------------------------

mat1 <- results1[[4]]$consensusMatrix
df1 <- data.frame(Class = Pdata1$`diagnosis:ch1`,
                  Consensus_class=as.character(results1[[4]]$consensusClass))

h1 <- HeatmapAnnotation(df = df1)
dimnames(mat1) <- list(colnames(Pdata1$`diagnosis:ch1`), colnames(Pdata1$`diagnosis:ch1`))
Heatmap(mat1,
        top_annotation = h1,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_title = "Consensus Clustering of GSE104948 most variable genes")



mat2 <- results2[[4]]$consensusMatrix
df2 <- data.frame(Class = Pdata2$disease,
                  Consensus_class=as.character(results2[[4]]$consensusClass))

h2 <- HeatmapAnnotation(df = df2)
dimnames(mat2) <- list(colnames(Pdata2$disease), colnames(Pdata2$disease))
Heatmap(mat2,
        top_annotation = h2,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_title = "Consensus Clustering of GSE37460 most variable genes")



mat3 <- results3[[4]]$consensusMatrix
df3 <- data.frame(Class = Pdata3$`disease status:ch1`,
                  Consensus_class=as.character(results3[[4]]$consensusClass))

h3 <- HeatmapAnnotation(df = df3)
dimnames(mat3) <- list(colnames(Pdata3$`disease status:ch1`), colnames(Pdata3$`disease status:ch1`))
Heatmap(mat3,
        top_annotation = h3,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_title = "Consensus Clustering of GSE32591 most variable genes")


####### PCA ############

###PCA using most variable genes------------------------------------
matdata1 <- as.matrix(data1)
pca1<- prcomp(t(matdata1), center=TRUE, scale. = TRUE)
ggbiplot(pca1,
         obs.scale = 1,
         var.scale = 1,
         var.axes = FALSE,
         groups = Pdata1$`diagnosis:ch1`) +
  ggtitle("PCA of GSE104948")



matdata2 <- as.matrix(data2)
pca2<- prcomp(t(matdata2), center=TRUE, scale. = TRUE)
ggbiplot(pca2,
         obs.scale = 1,
         var.scale = 1,
         var.axes = FALSE,
         groups = Pdata2$disease) +
  ggtitle("PCA of GSE37460")



matdata3 <- as.matrix(data3)
pca3<- prcomp(t(matdata3), center=TRUE, scale. = TRUE)
ggbiplot(pca3,
         obs.scale = 1,
         var.scale = 1,
         var.axes = FALSE,
         groups = Pdata3$disease) +
  ggtitle("PCA of GSE32591")



