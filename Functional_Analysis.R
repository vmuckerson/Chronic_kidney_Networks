##Load Libraries
library(msigdbr)
library(SPIA)
library(tidyverse)
library(dplyr)
library(plyr)
library(progeny)



## Load Data
setwd("/Users/saezlab/Documents/CKD_Data/Glom_Data/")

data1 <- readRDS(file = "./all_104948.rds")
Pdata1 <- readRDS(file = "./Pdata_104948.rds")
fdata1 <- readRDS(file = "./fdata_104948.rds")

data3 <- readRDS(file = "./all_GSE32591.rds")
Pdata3 <- readRDS(file = "./Pdata_GSE32591.rds")
fdata3 <- readRDS(file = "./fdata_GSE32591.rds")

DEA_results1 <- readRDS(file="./LN_DEA_104948.rds")
DEA_results3 <- readRDS(file="./LN_DEA_32591.rds")


## Select DE genes between LN and LD phenotypes

limma_genes <- filter(DEA_results1, adj.P.Val < 0.05)[[1]]
limma_genes3 <- filter(DEA_results3, adj.P.Val < 0.05)[[1]]

DEgenes1 <- DEA_results1$Gene == limma_genes
DEG1 <- DEA_results1[DEgenes1,]

DEgenes3 <- DEA_results3$Gene == limma_genes3
DEG3 <- DEA_results3[DEgenes3,]


## Annotate DE genes by Entrez ID instead of Gene symbol

matching_table1 <- fdata1[,c("ENTREZ_GENE_ID", "Symbol")] %>%
  as_tibble() %>%
  mutate(ENTREZ_GENE_ID = as.character(ENTREZ_GENE_ID)) %>%
  dplyr::rename(Gene = Symbol)

df1 <- DEG1 %>%
  rownames_to_column("Gene") %>%
  as_tibble()

joined_df <- inner_join(matching_table1, DEG1, by="Gene") %>%
  data.frame(row.names = 1, check.names = F, stringsAsFactors = F) %>%
  as.data.frame()

all1 <- inner_join(matching_table1, DEA_results1, by="Gene") %>%
  data.frame(row.names = 1, check.names = F, stringsAsFactors = F) %>%
  as.data.frame() %>% name_rows() %>% dplyr::select(-Gene, -logFC, -AveExpr, -t, -P.Value, -adj.P.Val, -B)


matching_table3 <- fdata3[,c("ID", "Gene Symbol")] %>%
  as_tibble() %>%
  mutate(ID = as.character(ID)) %>%
  dplyr::rename(Gene = "Gene Symbol")

df3 <- DEG3 %>%
  rownames_to_column("Gene") %>%
  as_tibble()

joined_df3 <- inner_join(matching_table3, DEG3, by="Gene") %>%
  data.frame(row.names = 1, check.names = F, stringsAsFactors = F) %>%
  as.data.frame()

all3 <- inner_join(matching_table3, DEA_results3, by="Gene") %>%
  data.frame(row.names = 1, check.names = F, stringsAsFactors = F) %>%
  as.data.frame() %>% name_rows() %>% dplyr::select(-Gene, -logFC, -AveExpr, -t, -P.Value, -adj.P.Val, -B)


de1 <- joined_df %>% dplyr::select(logFC) %>% t()
names(de1) <- colnames(de1)

de3 <- joined_df3 %>% dplyr::select(logFC) %>% t()
names(de3) <- colnames(de3)



## MSig DB to obtain reference gene set and pathway information

msig <- msigdbr(species = "Homo sapiens", category = NULL,
        subcategory = NULL)


## SPIA

spia_results1 <- spia(de = de1,
                      all = all1$.rownames,
                      organism = "hsa",
                      data.dir = NULL,
                      pathids = NULL,
                      nB = 2000,
                      plots = TRUE,
                      combine = "norminv")

spia_results3 <- spia(de = de3,
                      all = all3$.rownames,
                      organism = "hsa",
                      data.dir = NULL,
                      pathids = NULL,
                      nB = 2000,
                      plots = TRUE,
                      combine = "norminv")

## ORA performed via PANTHER to compare pathway predictions and choose the best one for further analysis. SPIA predicted SLE so I went with this one.

## Progeny

pathways <- progeny(as.matrix(data3), scale = FALSE)
controls <- Pdata3$disease == "LD"
ctmean <- apply(pathways[controls,], 2, mean)
ctsd <- apply(pathways[controls,], 2, sd)
ctmed <- apply(pathways[controls,], 2, median)
pathways_mean <- t(apply(pathways, 1, function(x) x - ctmean))
pathways_sd <- apply(pathways, 1, function(x) x/ctsd)
pathways_med <- apply(pathways, 1, function(x) x - ctmed)

Pln <- Pdata3$disease == "LN"
ln <- data3[Pln,]

samples <- intersect(colnames(data3), rownames(pathways))
sle <- data3[,samples]
mapk <- pathways[samples, "MAPK"]

association <- lm(ctmean ~ mapk)
