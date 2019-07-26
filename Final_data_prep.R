#Download libraries

library(GEOquery)
library(biomaRt)
library(qusage)
library(Biobase)
library(annotate)
library(hgu133a.db)
library(ggbiplot)

setwd("/Users/saezlab/Documents/CKD_Data/FinalData/GSE104948:54/")

###Download the data of the datasets as an Expression Set###########

Glom_data <- getGEO(GEO = 'GSE104948')
Glom_22945 <- Glom_data$`GSE104948-GPL22945_series_matrix.txt.gz`
Glom_24120 <- Glom_data$`GSE104948-GPL24120_series_matrix.txt.gz`

Tub_data <- getGEO(GEO = 'GSE104954')
Tub_22945 <- Tub_data$`GSE104954-GPL22945_series_matrix.txt.gz`
Tub_24120 <- Tub_data$`GSE104954-GPL24120_series_matrix.txt.gz`

#####Separate the phenotypic and expression data from the large files#####

#Phenotypic
PGlom_22945 <- pData(Glom_22945)
PGlom_24120 <- pData(Glom_24120)

PTub_22945 <- pData(Tub_22945)
PTub_24120 <- pData(Tub_24120)

#Expression data as a data frame
ExGlom_22945 <- as.data.frame(exprs(Glom_22945))
ExGlom_24120 <- as.data.frame(exprs(Glom_24120))
ExTub_22945 <- as.data.frame(exprs(Tub_22945))
ExTub_24120 <- as.data.frame(exprs(Tub_24120))


#####Annotate data frames with Gene Symbols instead of IDs#######

fGlom_22945 <- fData(Glom_22945)
fTub_22945 <- fData(Tub_22945)
fGlom24120 <- fData(Glom_24120)
fTub24120 <- fData(Tub_24120)

##Data without Gene Symbol data

genesymbols = getSYMBOL(x = as.character(fGlom24120$ENTREZ_GENE_ID), data = "hgu133a.db")


#Converting probe IDs to Gene Symbols using BiomaRT
mart1 <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

query <- fGlom24120$ENTREZ_GENE_ID
my_attributes <- c('hgnc_symbol', 'entrezgene_id')
matching_tableG24 <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                      filters = "entrezgene_id",
                      values = query,
                      bmHeader = T,
                      mart = mart1)


query <- fTub24120$ENTREZ_GENE_ID
matching_tableT24 <- getBM(my_attributes,
                          filters = 'entrezgene_id',
                          values = query,
                          mart = mart1)

colnames(matching_tableG24) <- c('probe', 'gene')
colnames(matching_tableT24) <- c('probe', 'gene')

matching_tableG24[matching_tableG24$gene == "NA"] <- NA
matching_tableT24[matching_tableT24$gene == "NA"] <- NA
is.na(matching_tableG24)
is.na(matching_tableT24)

#Data with Gene Symbol data
matching_tableG22 <- fGlom_22945[,c("ID", "Symbol")]
matching_tableT22 <- fTub_22945[,c("ID", "Symbol")]

dfG22 <- ExGlom_22945
dfG24 <- ExGlom_24120
dfT22 <- ExTub_22945
dfT24 <- ExTub_24120

dfG24$gene <- fGlom24120$ENTREZ_GENE_ID
rownames(dfG24) = as.character((dfG24$gene))

dfT24$gene <- fTub24120$ENTREZ_GENE_ID
rownames(dfT24) = as.character((dfT24$gene))

#duplicated entrez were discarded
dup_entrez = matching_tableG24$gene[duplicated(matching_tableG24$`HGNC symbol`)]
matching_tableG24 = matching_tableG24[!matching_tableG24$`HGNC symbol` %in% (dup_entrez),]
dfG24 = dfG24[as.character(matching_tableG24$`HGNC symbol`),]

dup_entrez2 = matching_tableT24$hgnc_symbol[duplicated(matching_tableT24$gene)]
matching_tableT24 = matching_tableT24[!matching_tableT24$hgnc_symbol %in% (dup_entrez2),]
dfT24 = dfT24[as.character(matching_tableT24$hgnc_symbol),]

#Adding Gene symbols to df
dfG22$gene <- matching_tableG22$Symbol
dfT22$gene <- matching_tableT22$Symbol
dfG24$gene <- matching_tableG24$`HGNC symbol`
dfT24$gene <- matching_tableT24$hgnc_symbol


#Aggregate the probe values
df2G22 <- aggregate(. ~gene, data = dfG22, mean)
df2T22 <- aggregate(. ~gene, data = dfT22, mean)
df2G24 <- aggregate(. ~gene, data = dfG24, mean)
df2T24 <- aggregate(. ~gene, data = dfT24, mean)

#annotate df
row.names(df2G22) <- dfG22$gene
row.names(df2T22) <- dfT22$gene
df2G22 <- df2G22[,-1]
df2T22 <- df2T22[,-1]
colnames(df2G22) <- Glom_22945@phenoData@data$geo_accession
colnames(df2T22) <- Tub_22945@phenoData@data$geo_accession

row.names(df2G24) <- df2G24$gene
row.names(df2T24) <- df2T24$gene
df2G24 <- df2G24[,-1]
df2T24 <- df2T24[,-1]
colnames(df2G24) <- Glom_24120@phenoData@data$geo_accession
colnames(df2T24) <- Tub_24120@phenoData@data$geo_accession


#Save everything
saveRDS(df2G22, file="./allGlom_22945.rds")
saveRDS(df2G24, file="./allGlom_24120.rds")
saveRDS(df2T22, file="./allTub_22945.rds")
saveRDS(df2G24, file="./allTub_24120.rds")

saveRDS(PGlom_22945, file="./PGlom_22945.rds")
saveRDS(PGlom_24120, file="./PGlom_24120.rds")
saveRDS(PTub_22945, file="./PTub_22945.rds")
saveRDS(PTub_24120, file="./PTub_24120.rds")

saveRDS(ExGlom_22945, file="./ExGlom_22945.rds")
saveRDS(ExGlom_24120, file="./ExGlom_24120.rds")
saveRDS(ExTub_22945, file="./ExTub_22945.rds")
saveRDS(ExTub_24120, file="./ExTub_24120.rds")






##################
##PCAs

#Glomerular data set
matGlom22 <- as.matrix(Glom_22945)
Gpca22<- prcomp(t(matGlom22), center=TRUE, scale. = TRUE)
ggbiplot(Gpca22,
         obs.scale = 1,
         var.scale = 1,
         var.axes = FALSE,
         groups = PGlom_22945$`diagnosis:ch1`) +
  ggtitle("PCA of GSE104948 Platform GPL22945")


matGlom24 <- as.matrix(Glom_24120)
Gpca24<- prcomp(t(matGlom24), center=TRUE, scale. = TRUE)
ggbiplot(Gpca24,
         obs.scale = 1,
         var.scale = 1,
         var.axes = FALSE,
         groups = PGlom_24120$`diagnosis:ch1`) +
  ggtitle("PCA of GSE104948 Platform GPL24120")


#Tubular Data Set

matTub22 <- as.matrix(Tub_22945)
Tpca22<- prcomp(t(matTub22), center=TRUE, scale. = TRUE)
ggbiplot(Tpca22,
         obs.scale = 1,
         var.scale = 1,
         var.axes = FALSE,
         groups = PTub_22945$`diagnosis:ch1`) +
  ggtitle("PCA of GSE104954 Platform GPL22945")

matTub24 <- as.matrix(Tub_24120)
Tpca24 <- prcomp(t(matTub24), center=TRUE, scale. = TRUE)
ggbiplot(Tpca24,
         obs.scale = 1,
         var.scale = 1,
         var.axes = FALSE,
         groups = PTub_24120$`diagnosis:ch1`) +
  ggtitle("PCA of GSE10954 Platform GPL24120")



saveRDS(Gpca22, file="./Gpca22.rds")
saveRDS(Gpca24, file="./Gpca24.rds")
saveRDS(Tpca22, file="./Tpca22.rds")
saveRDS(Tpca24, file="./Tpca24.rds")


library(Rtsne)


#Multidimensional scaling (mds)

gsne22 <- Rtsne(t(matTub22), perplexity = 10, theta = 0.01)

plot(gsne22$Y, col = as.factor(Tub_22945$`diagnosis:ch1`), pch = 19,
     main = "tSNE plot for GSE104954 Platform 22945")










