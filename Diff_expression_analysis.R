library(limma)
library(tidyverse)
library(BiocManager)

##### Diff Expression Analysis LN to Healthy (GEO32591) ##########

#Download data

my_data <- readRDS(file = "./all_GSE32591.rds")
Pdata <- readRDS(file = "./Pdata_GSE32591.rds")


#Subset the data by condition (Lupus Nephritis or Living Donor)

f = factor(Pdata$disease, levels = c("LD", "LN"))
design <- model.matrix(~0+f)
colnames(design) <- c("LD", "LN")


#Differential Expression Analysis with Limma

fit <- lmFit(my_data, design)
cont.matrix <- makeContrasts(Condition_dif = LN - LD, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)


#Selecting significant genes from Limma results

DEA_results <- as.data.frame(topTable(fit2,
                                      adjust.method = "BH",
                                      number = Inf)) %>% rownames_to_column(var = "Gene") %>%
               arrange(desc(abs(t))) %>% as_tibble()

limma_genes <- filter(DEA_results, adj.P.Val < 0.05)[[1]]






##### Diff Expression Analysis LN to Healthy from 2nd dataset (GEO104948) ##########

#Download data

my_data2 <- readRDS(file = "./all_104948.rds")
Pdata2 <- readRDS(file = "./Pdata_104948.rds")


#Subset the data by condition (Lupus Nephritis or Living Donor)
Pdata2$`diagnosis:ch1`[is.na(Pdata2$`diagnosis:ch1`)] = "Healthy"

f2 <- factor(Pdata2$`diagnosis:ch1`, levels = c("Systemic Lupus Erythematosus", "Healthy"))
design2 <- model.matrix(~0+f2)
colnames(design2) <- c("SLE", "Healthy")


#Subsetting data for only SLE and Healthy samples
Pdata2$disease <- Pdata2$`diagnosis:ch1`
mine <- (Pdata2$disease == "Systemic Lupus Erythematosus" | Pdata2$disease == "Healthy")
limma_data2 <- my_data2[,mine]


#Differential Expression Analysis with Limma between SLE and Healthy samples
fit3 <- lmFit(limma_data2, design2)
cont.matrix2 <- makeContrasts(Condition_dif = SLE - Healthy, levels = design2)
fit4 <- contrasts.fit(fit3, cont.matrix2)
fit4 <- eBayes(fit3)


#Selecting significant genes from Limma results
DEA_results2 <- as.data.frame(topTable(fit4,
                                      adjust.method = "BH",
                                      number = Inf)) %>% rownames_to_column(var = "Gene") %>%
  arrange(desc(abs(F))) %>% as_tibble()

limma_genes2 <- filter(DEA_results2, adj.P.Val < 0.05)[[1]]





#Selecting genes that are sig dif exp in both DEAs

footprint <- limma_genes2 %in% limma_genes
footprint <- limma_genes2[footprint]








