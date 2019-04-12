library(tidyverse)
library(ggplot2)
library(GEOquery)

A_N570 <- readRDS("/Users/saezlab/Documents/CKD_Data/Norm_Data/GSE69438_Norm570.rds")
Apheno <- getGEO(GEO = "GSE69438", destdir = tempdir("/Users/saezlab/Documents/CKD_Data/GSE69438/"))

ADisease <- Apheno$GSE69438_series_matrix.txt.gz$title %>% as.tibble %>%
  separate(value, into = c("Disease", "Combo"), sep = " Cprobe") %>%
  select(-c(Combo)) %>% print()

ATissue <- Apheno$GSE69438_series_matrix.txt.gz$characteristics_ch1 %>% as.tibble %>%
  separate(value, into = c("t", "Tissue", "y", "z", "j"), sep = " ") %>%
  select(-c(t, y, z, j)) %>%
  print()

#AGene <- Apheno$GSE69438_series_matrix.txt.gz@featureData@data$Description %>% as.tibble() %>% print()

v <- tibble(sample = colnames(A_N570)[-1], tissue = ATissue$Tissue, disease = ADisease$Disease)
v




B_N570 <- readRDS("/Users/saezlab/Documents/CKD_Data/Norm_Data/GSE104948_Norm570.rds")
Bpheno <- getGEO(GEO = "GSE104948", destdir = tempdir("/Users/saezlab/Documents/CKD_Data/GSE104948/"))

B <- Bpheno$`GSE104948-GPL22945_series_matrix.txt.gz`$`title`%>%
  as.tibble() %>%
  separate(value, into = c("H", "tissue", "disease"), sep = "-") %>%
  select(-c(H)) %>%
  mutate(disease= str_extract(disease, '\\D*(?=\\d)')) %>%
           print()

#BGene <- Bpheno$`GSE104948-GPL22945_series_matrix.txt.gz`@featureData@data$SPOT_ID %>% as.tibble() %>% print()

w <- tibble(sample = colnames(B_N570)[-1], tissue = B$tissue, disease = B$disease)
w




B_N96 <- readRDS("/Users/saezlab/Documents/CKD_Data/Norm_Data/GSE104948_Norm96.rds")

B <- Bpheno$`GSE104948-GPL24120_series_matrix.txt.gz`$title %>% as.tibble %>%
  separate(value, into = c("h", "tissue", "Disease"), sep = "-") %>%
  mutate(disease= str_extract(Disease, '\\D*(?=\\d)')) %>%
  select(-c(h, Disease)) %>%
  print()

#BGene96 <- Bpheno$`GSE104948-GPL24120_series_matrix.txt.gz`@featureData@data$Description %>% as.tibble() %>% print()

BSample96 <- Bpheno$`GSE104948-GPL24120_series_matrix.txt.gz`@phenoData@data$geo_accession
BSample96

x <- tibble(sample = BSample96[1:125], tissue = B$tissue, disease = B$disease)
x




C_N570 <- readRDS("/Users/saezlab/Documents/CKD_Data/Norm_Data/GSE104954_Norm570.rds")
Cpheno <- getGEO(GEO = "GSE104954", destdir = tempdir("/Users/saezlab/Documents/CKD_Data/GSE104954/"))

C <- Cpheno$`GSE104954-GPL22945_series_matrix.txt.gz`$title %>% as.tibble %>%
  separate(value, into = c("h", "tissue", "Disease"), sep = "-") %>%
  mutate(disease= str_extract(Disease, '\\D*(?=\\d)')) %>%
  select(-c(h, Disease)) %>%
  print()

#CGene <- Cpheno$`GSE104954-GPL22945_series_matrix.txt.gz`@featureData@data$SPOT_ID %>% as.tibble() %>% print()

y <- tibble(sample = colnames(C_N570)[-1], tissue = C$tissue, disease = C$disease)
y




C_N96 <- readRDS("/Users/saezlab/Documents/CKD_Data/Norm_Data/GSE104954_Norm96.rds")

C <- Cpheno$`GSE104954-GPL24120_series_matrix.txt.gz`$title %>% as.tibble %>%
  separate(value, into = c("h", "tissue", "Disease"), sep = "-") %>%
  mutate(disease= str_extract(Disease, '\\D*(?=\\d)')) %>%
  select(-c(h, Disease)) %>%
  print()

#CGene96 <- Cpheno$`GSE104954-GPL24120_series_matrix.txt.gz`@featureData@data$SPOT_ID %>% as.tibble() %>% print()

z <- tibble(sample = colnames(C_N96)[-1], tissue = C$tissue, disease = C$disease)
z




