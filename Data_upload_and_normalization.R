#Download libraries

library(GEOquery)
library(biomaRt)
library(qusage)
library(Biobase)
library(annotate)
library(hgu133a.db)

setwd("/Users/saezlab/Documents/CKD_Data/FinalData/")

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
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

query <- fGlom24120$ENTREZ_GENE_ID
my_attributes <- c('hgnc_symbol', 'entrezgene_id')
matching_tableG24 <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                      filters = "entrezgene_id",
                      values = query,
                      bmHeader = T,
                      mart = mart)


query <- fTub24120$ENTREZ_GENE_ID
matching_tableT24 <- getBM(my_attributes,
                          filters = 'entrezgene_id',
                          values = query,
                          mart = mart)

colnames(matching_tableG24) <- c('probe', 'gene')
colnames(matching_tableT24) <- c('probe', 'gene')

matching_tableG24$gene[matching_tableG24$gene == "NA"] <- NA
matching_tableT24$gene[matching_tableT24$gene == "NA"] <- NA
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
matching_tableG24 = matching_tableG24[matching_tableG24$probe!="",]

dfT24$gene <- fTub24120$ENTREZ_GENE_ID
rownames(dfT24) = as.character((dfT24$gene))
matching_tableT24 = matching_tableT24[matching_tableT24$probe!="",]

#duplicated entrez were discarded
dup_entrez = matching_tableG24$gene[duplicated(matching_tableG24$gene)]
matching_tableG24 = matching_tableG24[!matching_tableG24$gene %in% (dup_entrez),]
dfG24 = dfG24[as.character(matching_tableG24$gene),]

dup_entrez2 = matching_tableT24$gene[duplicated(matching_tableT24$gene)]
matching_tableT24 = matching_tableT24[!matching_tableT24$gene %in% (dup_entrez2),]
dfT24 = dfT24[as.character(matching_tableT24$gene),]

#Adding Gene symbols to df
dfG22$gene <- matching_tableG22$Symbol
dfT22$gene <- matching_tableT22$Symbol
dfG24$gene <- matching_tableG24$probe
dfT24$gene <- matching_tableT24$probe


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
saveRDS(df2G22, file="./GSE104948/allGlom_22945.rds")
saveRDS(df2G24, file="./GSE104948/allGlom_24120.rds")
saveRDS(df2T22, file="./GSE104948/allTub_22945.rds")
saveRDS(df2G24, file="./GSE104948/allTub_24120.rds")

saveRDS(PGlom_22945, file="./GSE104948/PGlom_22945.rds")
saveRDS(PGlom_24120, file="./GSE104948/PGlom_24120.rds")
saveRDS(PTub_22945, file="./GSE104948/PTub_22945.rds")
saveRDS(PTub_24120, file="./GSE104948/PTub_24120.rds")

saveRDS(ExGlom_22945, file="./GSE104948/ExGlom_22945.rds")
saveRDS(ExGlom_24120, file="./GSE104948/ExGlom_24120.rds")
saveRDS(ExTub_22945, file="./GSE104948/ExTub_22945.rds")
saveRDS(ExTub_24120, file="./GSE104948/ExTub_24120.rds")

