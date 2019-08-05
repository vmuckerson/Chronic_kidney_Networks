library(limma)
library(tidyverse)
library(BiocManager)
library(ConsensusClusterPlus)
library(ComplexHeatmap)

##### Diff Expression Analysis LN to Healthy (GEO32591) ##########

#Download data
setwd("/Users/saezlab/Documents/CKD_Data/Glom_Data/")
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
fit4 <- eBayes(fit4)


#Selecting significant genes from Limma results
DEA_results2 <- as.data.frame(topTable(fit4,
                                      adjust.method = "BH",
                                      number = Inf)) %>% rownames_to_column(var = "Gene") %>%
  arrange(desc(abs(t))) %>% as_tibble()
limma_genes2 <- filter(DEA_results2, adj.P.Val < 0.05)[[1]]


#Selecting genes that are sig dif exp in both DEAs

footprint_genes <- limma_genes2 %in% limma_genes
footprint <- limma_genes2[footprint_genes]






############# Matrix Scoring #################

#limiting t values to those present in DEGs of interest
dea2 <- DEA_results2$Gene %in% limma_genes2
matB <- DEA_results2[dea2,]
matB <- na.omit(matB)

#limiting 2nd matrix t values to genes matching those in 1st matrix
dea <- DEA_results$Gene %in% limma_genes2
matA <- DEA_results[dea,]
matA <- na.omit(matA)

#limiting data of 2nd matrix to the genes of interest
limit <- rownames(limma_data2) %in% limma_genes2
risk_data2 <- limma_data2[limit,]
risk_data2 <- na.omit(risk_data2)

#limiting data of 1st matrix to the genes of interest
edit <- rownames(my_data) %in% matB$Gene
risk_data <- my_data[edit,]
risk_data <- risk_data[matB$Gene,]
risk_data <- na.omit(risk_data)

#limiting matrices again to ensure the same size
fix <- matA$Gene %in% rownames(risk_data2)
matA<- matA[fix,]

fix2 <- matB$Gene %in% rownames(risk_data)
matB<- matB[fix2,]

#matrix multiplication of one mat with the other's t values
risk <- matB$t %*% as.matrix(risk_data[matB$Gene,])


#limiting pdata of 2nd matrix to the samples chosen of interest
rip <- Pdata$geo_accession %in% colnames(risk_data)
ripdata <- Pdata[rip,]
ripdata <- na.omit(ripdata)

rip2 <- Pdata2$geo_accession %in% colnames(risk_data2)
ripdata2 <- Pdata2[rip2,]
ripdata2 <- na.omit(ripdata2)

risk <- t(risk)

riskmat <- as.tibble(risk) %>%
  add_column(as.factor(ripdata$disease)) %>% as.data.frame()
colnames(riskmat) = c("Value", "Disease")

#creating a boxplot of healthy vs sick values
boxplot(data = riskmat,
        Value~Disease,
        main = "Boxplot of Healthy vs. LN Risk Index for GSE32591 using DEGs from GSE104948",
        cex.main = 0.85,
        col = c("turquoise", "orange"))
#legend(0.75, -7500, c("Healthy", "LN"), fill = c("turquoise", "orange"))

