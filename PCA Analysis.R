library(tidyverse)
library(ggplot2)
library(GEOquery)
library(ggbiplot)
library(BiocManager)



A_N570 <- readRDS("/Users/saezlab/Documents/CKD_Data/Norm_Data/GSE69438_Norm570.rds")
rownames(A_N570) <- A_N570$GeneSymbol
Apheno <- getGEO(GEO = "GSE69438", destdir = tempdir("/Users/saezlab/Documents/CKD_Data/GSE69438/"))

Adisease <- Apheno$GSE69438_series_matrix.txt.gz$title %>% as.tibble %>%
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

Atissue <- Apheno$GSE69438_series_matrix.txt.gz$characteristics_ch1 %>% as.tibble %>%
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
  add_column(disease = Adisease, .before = 1) %>%
  add_column(tissue = Atissue, .before = 1) %>% 
  add_column("experiment" = "GSE69438", .before = 1) %>%
  add_column("platform" = "GPL570", .after = 1)


A.pca<- prcomp(AEx[,c(5:23525)], center=TRUE, scale. = TRUE)
ggbiplot(A.pca, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups=Adisease) +
  ggtitle("PCA of GSE69438 Platform GPL570")



B_N570 <- readRDS("/Users/saezlab/Documents/CKD_Data/Norm_Data/GSE104948_Norm570.rds")
rownames(B_N570) <- B_N570$GeneSymbol
Bpheno <- getGEO(GEO = "GSE104948", destdir = tempdir("/Users/saezlab/Documents/CKD_Data/GSE104948/"))

B <- Bpheno$`GSE104948-GPL22945_series_matrix.txt.gz`$`title`%>%
  as.tibble() %>%
  separate(value, into = c("H", "tissue", "disease"), sep = "-") %>%
  select(-c(H)) %>%
  mutate(disease= str_extract(disease, '[a-zA-Z]+')) %>%
           print()

BEx <- t(B_N570[-1]) %>% as.tibble %>%
  add_column(disease = B$disease, .before = 1) %>%
  add_column(tissue = B$tissue, .before = 1) %>%
  add_column("experiment" = "GSE104948", .before = 1)%>%
  add_column("platform" = "GPL570", .after = 1)

B.pca570<- prcomp(BEx[,c(5:23525)], center=TRUE, scale. = TRUE)
ggbiplot(B.pca570, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups=B$disease) +
  ggtitle("PCA of GSE104948 Platform GPL570")



B_N96 <- readRDS("/Users/saezlab/Documents/CKD_Data/Norm_Data/GSE104948_Norm96.rds")
rownames(B_N96) <- B_N96$Group.1
B96 <- Bpheno$`GSE104948-GPL24120_series_matrix.txt.gz`$title %>% as.tibble %>%
  separate(value, into = c("h", "tissue", "Disease"), sep = "-") %>%
  mutate(disease= str_extract(Disease, '[a-zA-Z]+')) %>%
  select(-c(h, Disease)) %>%
  print()

BEx_96 <- t(B_N96[-1]) %>% as.tibble %>%
  add_column(disease = B96$disease, .before = 1) %>%
  add_column(tissue = B96$tissue, .before = 1) %>%
  add_column("experiment" = "GSE104948", .before = 1)%>%
  add_column("platform" = "GPL96", .after = 1)

B.pca96<- prcomp(BEx_96[,c(5:13520)], center=TRUE, scale. = TRUE)
ggbiplot(B.pca96, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups=B96$disease) +
  ggtitle("PCA of GSE104948 Platform GPL96")



C_N570 <- readRDS("/Users/saezlab/Documents/CKD_Data/Norm_Data/GSE104954_Norm570.rds")
rownames(C_N570) <- C_N570$GeneSymbol
Cpheno <- getGEO(GEO = "GSE104954", destdir = tempdir("/Users/saezlab/Documents/CKD_Data/GSE104954/"))

C <- Cpheno$`GSE104954-GPL22945_series_matrix.txt.gz`$title %>% as.tibble %>%
  separate(value, into = c("h", "tissue", "Disease"), sep = "-") %>%
  mutate(disease= str_extract(Disease, '[a-zA-Z]+')) %>%
  select(-c(h, Disease)) %>%
  print()

CEx <- t(C_N570[-1]) %>% as.tibble %>%
  add_column(disease = C$disease, .before = 1) %>%
  add_column(tissue = C$tissue, .before = 1) %>%
  add_column("experiment" = "GSE104954", .before = 1)%>%
  add_column("platform" = "GPL570", .after = 1)

C.pca570<- prcomp(CEx[,c(5:23525)], center=TRUE, scale. = TRUE)
ggbiplot(C.pca570, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups=C$disease) +
  ggtitle("PCA of GSE104954 Platform GPL570")



C_N96 <- readRDS("/Users/saezlab/Documents/CKD_Data/Norm_Data/GSE104954_Norm96.rds")
rownames(C_N96) <- C_N96$GeneSymbol
C <- Cpheno$`GSE104954-GPL24120_series_matrix.txt.gz`$title %>% as.tibble %>%
  separate(value, into = c("h", "tissue", "Disease"), sep = "-") %>%
  mutate(disease= str_extract(Disease, '[a-zA-Z]+')) %>%
  select(-c(h, Disease)) %>%
  print()

CEx_96 <- t(C_N96[-1]) %>% as.tibble %>%
  add_column(disease = C$disease, .before = 1) %>%
  add_column(tissue = C$tissue, .before = 1) %>%
  add_column("experiment" = "GSE104954", .before = 1)%>%
  add_column("platform" = "GPL570", .after = 1)

C.pca96<- prcomp(CEx_96[,c(5:13520)], center=TRUE, scale. = TRUE)
ggbiplot(C.pca96, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups=C$disease) +
  ggtitle("PCA of GSE10954 Platform GPL96")


AEx %>% as.list()
BEx %>% as.list()
CEx %>% as.list()
list_570 <- list(AEx, BEx, CEx)
idx <- Reduce(intersect, lapply(list_570, colnames)) %>% as_vector()

A570 <- AEx[, idx]
B570 <- BEx[, idx]
C570 <- CEx[, idx]

all570 <- rbind(A570, B570, C570)
all570.pca<- prcomp(all570[,c(5:173)], center=TRUE, scale. = TRUE)

ggbiplot(all570.pca, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups=all570$experiment) +
  ggtitle("PCA of Platform 570") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(colour = "Experiment")

ggbiplot(all570.pca, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups=all570$tissue) +
  ggtitle("PCA of Platform 570") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(colour = "Tissue Type")

ggbiplot(all570.pca, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups=all570$disease) +
  ggtitle("PCA of Platform 570") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(colour = "Disease")



BEx_96 %>% as.list()
CEx_96 %>% as.list()
list_96 <- list(BEx_96, CEx_96)
idx_96 <- Reduce(intersect, lapply(list_96, colnames)) %>% as_vector()

B96 <- BEx_96[, idx_96]
C96 <- CEx_96[, idx_96]

all96 <- rbind(B96, C96)
all96.pca<- prcomp(all96[,c(5:260)], center=TRUE, scale. = TRUE)

ggbiplot(all96.pca, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups=all96$experiment) +
  ggtitle("PCA of Platform 96") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(colour = "Experiment")

ggbiplot(all96.pca, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups=all96$tissue) +
  ggtitle("PCA of Platform 96") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(colour = "Tissue Type")

ggbiplot(all96.pca, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups=all96$disease) +
  ggtitle("PCA of Platform 96") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(colour = "Disease")


list_all <- list(AEx, BEx, BEx_96, CEx, CEx_96)
idx_all <- Reduce(intersect, lapply(list_all, colnames)) %>% as_vector()

Ax <- AEx[, idx_all]
Bx <- BEx[, idx_all]
By <- BEx_96[, idx_all]
Cx <- CEx[, idx_all]
Cy <- CEx_96[, idx_all]

all_data <- rbind(Ax, Bx, By, Cx, Cy)
all.pca<- prcomp(all_data[,c(5:433)], center=TRUE, scale. = TRUE)

ggbiplot(all.pca, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups=all_data$experiment) +
  ggtitle("PCA of All Samples") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(colour = "Experiment")

ggbiplot(all.pca, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups=all_data$tissue) +
  ggtitle("PCA of All Samples") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(colour = "Tissue Type")

ggbiplot(all.pca, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups=all_data$disease) +
  ggtitle("PCA of All Samples") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(colour = "Disease")

ggbiplot(all.pca, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups=all_data$platform) +
  ggtitle("PCA of All Samples") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(colour = "Platform")

##################

all_data <- as.matrix(all_data)

batch_f <- all_data[,1] %>% as.factor()

comdata <- ComBat(all_data[c(5:433)], batch_f)

