#Download libraries

library(GEOquery)
library(biomaRt)
library(qusage)
library(Biobase)
library(annotate)
library(hgu133a.db)
library(ggbiplot)
library(tidyverse)

setwd("/Users/saezlab/Documents/CKD_Data/Glom_data/")

###Download the data of the datasets as an Expression Set###########

mydata <- getGEO(GEO = 'GSE104948')
data1 <- mydata[[1]]
data2 <- mydata2[[2]]


#####Separate the phenotypic, expression, and feature data from the large files#####

#Phenotypic
Pdata1 <- pData(data1)
Pdata2 <- pData(data2)


#Expression data as a data frame
Exdata1 <- as.data.frame(exprs(data1))
Exdata2 <- as.data.frame(exprs(data2))


#Feature data
fdata1 <- fData(data1)
fdata2 <- fData(data2)



#####Annotate data frames with Gene Symbols instead of IDs#######
##Data with Gene Symbol data
matching_table1 <- fdata1[,c("ID", "Symbol")] %>%
  as_tibble()

df1 <- Exdata1 %>%
  rownames_to_column("ID") %>%
  gather(sample, expression,-ID) %>%
  as_tibble()

df1 <- inner_join(matching_table1, df1, by="ID") %>%
  select(-ID) %>% 
  spread(sample, expression) %>%
  data.frame(row.names = 1, check.names = F, stringsAsFactors = F) %>%
  as.data.frame()

matching_table2 <- fdata1[,c("ID", "Symbol")] %>%
  as_tibble()

df2 <- Exdata2 %>%
  rownames_to_column("ID") %>%
  gather(sample, expression,-ID) %>%
  as_tibble()

df2 <- inner_join(matching_table2, df2, by="ID") %>%
  select(-ID) %>% 
  spread(sample, expression) %>%
  data.frame(row.names = 1, check.names = F, stringsAsFactors = F) %>%
  as.data.frame()


#Adding Gene symbols to df
df1$gene <- matching_table1$Symbol
df2$gene <- matching_table2$Symbol


#Aggregate the probe values
df12 <- aggregate(. ~gene, data = df1, mean)
df22 <- aggregate(. ~gene, data = df2, mean)

#annotate df
row.names(df12) <- df12$gene
df12 <- df12[,-1]
colnames(df12) <- Pdata1$geo_accession

row.names(df22) <- df22$gene
df22 <- df22[,-1]
colnames(df22) <- Pdata2$geo_accession

df12 <- na.omit(df12)
df22 <- na.omit(df22)


#Save everything
#saveRDS(df12, file="./all_104948.rds")
#saveRDS(Pdata1, file="./Pdata_104948.rds")
#saveRDS(Exdata1, file="./Ex_104948.rds")

#saveRDS(df22, file="./all_104954.rds")
#saveRDS(Pdata2, file="./Pdata_104954.rds")
#saveRDS(Exdata2, file="./Ex_104954.rds")




##################
##PCAs

#Glomerular data set
matdata1 <- as.matrix(data1)
pca1<- prcomp(t(matdata1), center=TRUE, scale. = TRUE)
ggbiplot(pca1,
         obs.scale = 1,
         var.scale = 1,
         var.axes = FALSE,
         groups = Pdata1$`diagnosis:ch1`) +
  ggtitle("PCA of GSE104948 Platform GPL22945")


matdata2 <- as.matrix(data2)
pca2<- prcomp(t(matdata2), center=TRUE, scale. = TRUE)
ggbiplot(pca2,
         obs.scale = 1,
         var.scale = 1,
         var.axes = FALSE,
         groups = Pdata2$`diagnosis:ch1`) +
  ggtitle("PCA of GSE104948 Platform GPL24120")

#saveRDS(pca1, file="./pca_104948.rds")
#saveRDS(pca2, file="./pca_104954.rds")


