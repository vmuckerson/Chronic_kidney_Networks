library(tidyverse)
library(ggplot2)
library(GEOquery)

A_N570 <- readRDS("/Users/saezlab/Documents/CKD_Data/Norm_Data/GSE69438_Norm570.rds")
Apheno <- getGEO(GEO = "GSE69438", destdir = tempdir("/Users/saezlab/Documents/CKD_Data/GSE69438/"))

disease <- Apheno$GSE69438_series_matrix.txt.gz$title %>% as.tibble %>%
  separate(value, into = c("Disease", "Combo"), sep = " Cprobe") %>%
  select(-c(Combo)) %>% as_vector %>% sapply(function(x) {
    x <- as.character(x)
    if(grepl("^Diabetic Nephropathy",x)) {
      "DN"
    } else if (grepl("^Focal and Segmental Glomerulosclerosis",x)) {
      "FSGS"
    } else if (grepl("^Hypertensive nephropathy",x)) {
      "HT"
    } else if (grepl("^IgA Nephropathy",x)) {
      "IgA"
    } else if (grepl("^Minimal Change Disease",x)) {
      "MCD"
    } else if (grepl("^Membranous Glomerulonephritis",x)) {
      "MGN"
    } else if (grepl("^CKD",x)) {
      "CKD"
    } else if (grepl("Lupus Nephritis",x)) {
      "LN"
    } else {
      NA
    }
  })

tissue <- Apheno$GSE69438_series_matrix.txt.gz$characteristics_ch1 %>% as.tibble %>%
  separate(value, into = c("t", "Tissue", "y", "z", "j"), sep = " ") %>%
  select(-c(t, y, z, j)) %>% as_vector() %>% sapply(function(x) {
    x <- as.character(x)
    if(grepl("^Tubulointerstitium",x)) {
      "Tub"
    } else {
      NA
    }
  })

AEx <- t(A_N570[-1]) %>% as.tibble %>%
  add_column(disease, .before = 1) %>%
  add_column(tissue, .before = 1) %>% 
  add_column("experiment" = "GSE69438", .before = 1) %>%
  add_column("platform" = "GPL570", .after = 1)

A.pca<- prcomp(AEx[,c(5:23525)], center=TRUE, scale. = TRUE)
ggbiplot(A.pca, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups=tissue) +
  ggtitle("PCA of GSE69438 Platform GPL570")



B_N570 <- readRDS("/Users/saezlab/Documents/CKD_Data/Norm_Data/GSE104948_Norm570.rds")
Bpheno <- getGEO(GEO = "GSE104948", destdir = tempdir("/Users/saezlab/Documents/CKD_Data/GSE104948/"))

B <- Bpheno$`GSE104948-GPL22945_series_matrix.txt.gz`$`title`%>%
  as.tibble() %>%
  separate(value, into = c("H", "tissue", "disease"), sep = "-") %>%
  select(-c(H)) %>%
  mutate(disease= str_extract(disease, '\\D*(?=\\d)')) %>%
           print()

BEx <- t(B_N570[-1]) %>% as.tibble %>%
  add_column("disease", B$disease, .before = 1) %>%
  add_column("tissue", B$tissue, .before = 1) %>%
  add_column("experiment" = "GSE104948", .before = 1)%>%
  add_column("platform" = "GPL570", .after = 1)

B.pca570<- prcomp(BEx[,c(5:23525)], center=TRUE, scale. = TRUE)
ggbiplot(B.pca570, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups=B$disease) +
  ggtitle("PCA of GSE104948 Platform GPL570")



B_N96 <- readRDS("/Users/saezlab/Documents/CKD_Data/Norm_Data/GSE104948_Norm96.rds")

B <- Bpheno$`GSE104948-GPL24120_series_matrix.txt.gz`$title %>% as.tibble %>%
  separate(value, into = c("h", "tissue", "Disease"), sep = "-") %>%
  mutate(disease= str_extract(Disease, '\\D*(?=\\d)')) %>%
  select(-c(h, Disease)) %>%
  print()

BEx_96 <- t(B_N96[-1]) %>% as.tibble %>%
  add_column("disease", B$disease, .before = 1) %>%
  add_column("tissue", B$tissue, .before = 1) %>%
  add_column("experiment" = "GSE104948", .before = 1)%>%
  add_column("platform" = "GPL96", .after = 1)

B.pca96<- prcomp(BEx_96[,c(5:13520)], center=TRUE, scale. = TRUE)
ggbiplot(B.pca96, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups=B$disease) +
  ggtitle("PCA of GSE104948 Platform GPL96")



C_N570 <- readRDS("/Users/saezlab/Documents/CKD_Data/Norm_Data/GSE104954_Norm570.rds")
Cpheno <- getGEO(GEO = "GSE104954", destdir = tempdir("/Users/saezlab/Documents/CKD_Data/GSE104954/"))

C <- Cpheno$`GSE104954-GPL22945_series_matrix.txt.gz`$title %>% as.tibble %>%
  separate(value, into = c("h", "tissue", "Disease"), sep = "-") %>%
  mutate(disease= str_extract(Disease, '\\D*(?=\\d)')) %>%
  select(-c(h, Disease)) %>%
  print()

CEx <- t(C_N570[-1]) %>% as.tibble %>%
  add_column("disease", C$disease, .before = 1) %>%
  add_column("tissue", C$tissue, .before = 1) %>%
  add_column("experiment" = "GSE104954", .before = 1)%>%
  add_column("platform" = "GPL570", .after = 1)

C.pca570<- prcomp(CEx[,c(5:23525)], center=TRUE, scale. = TRUE)
ggbiplot(C.pca570, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups=C$disease) +
  ggtitle("PCA of GSE104954 Platform GPL570")



C_N96 <- readRDS("/Users/saezlab/Documents/CKD_Data/Norm_Data/GSE104954_Norm96.rds")

C <- Cpheno$`GSE104954-GPL24120_series_matrix.txt.gz`$title %>% as.tibble %>%
  separate(value, into = c("h", "tissue", "Disease"), sep = "-") %>%
  mutate(disease= str_extract(Disease, '\\D*(?=\\d)')) %>%
  select(-c(h, Disease)) %>%
  print()

CEx_96 <- t(C_N96[-1]) %>% as.tibble %>% add_column(C$tissue, .before = 1) %>%
  add_column(C$disease, .before = 1) %>%
  add_column("experiment" = "GSE104954", .before = 1)%>%
  add_column("platform" = "GPL570", .after = 1)

C.pca96<- prcomp(CEx_96[,c(5:13520)], center=TRUE, scale. = TRUE)
ggbiplot(C.pca96, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups=C$disease) +
  ggtitle("PCA of GSE10954 Platform GPL96")

#alldata <- cbind(AEx, BEx, BEx_96, CEx, CEx_96)
#all.pca <- prcomp(alldata,center=TRUE, scale. = TRUE)
#ggbiplot(all.pca, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups=C$disease) +
#  ggtitle("PCA of All Datasets")
