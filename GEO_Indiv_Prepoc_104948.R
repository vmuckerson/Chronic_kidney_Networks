#Download libraries

library(GEOquery)
library(affy)
library(oligo)
library(Biobase)

wd <- file.path(getwd())
exper <- 'GSE104948'


### 1 Download the data in your local folder ###########
# Get the TAR file with all CEL files (each one sample)

#Download CEL files
url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104948/suppl/GSE104948_RAW.tar"
tar_file <- paste0("./CKD_Data/",basename(url))

# Create dir folder
dir.create("./CKD_Data")
download.file(url=url, destfile = tar_file)

# Untar files
untar(tar_file, exdir="./CKD_Data")


# get Metadata #######################
# Source (FTP): ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104948/matrix/

# 1st platform
url1 <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104948/matrix/GSE104948-GPL22945_series_matrix.txt.gz"
matrix_file1 <- paste0("./CKD_data/",basename(url1))
download.file(url1, matrix_file1)

# 2nd platform
url2 <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104948/matrix/GSE104948-GPL24120_series_matrix.txt.gz"
matrix_file2 <- paste0("./CKD_data/",basename(url2))
download.file(url2, matrix_file2)



### 2 Load the data ######################################
# Get all the CEL files
fls <- list.files(path = "./CKD_Data",pattern = "\\.CEL(\\.gz)$", full.names = TRUE)
names(fls) <- basename(fls)

# This dataset is a bit tricky. There are two platforms within the same dataset.
# Thus we cannot just load all the CEL files because they have different dimensions.
# Instead we have to load all samples by platform.



###################################################### GPL22945 ###################
GPL22945 <- getGEO(filename = "./CKD_Data/GSE104948-GPL22945_series_matrix.txt.gz")
GEOACC_GPL22945 <- basename(as.character(pData(GPL22945)$supplementary_file))

# Read raw data
rawData_GPL22945 <- oligo::read.celfiles(fls[GEOACC_GPL22945])
# Rename samples w/ only the GEO accession number
# Note: Not always this is needed
sampleNames(rawData_GPL22945) <- sapply(sampleNames(rawData_GPL22945),function(sampleName) strsplit(sampleName,split="_")[[1]][1])

# RMA normalization
norm_GPL22945 <- oligo::rma(rawData_GPL22945,background=TRUE, normalize=TRUE)

# Update sample metadata
stopifnot(all(sampleNames(norm_GPL22945) %in% sampleNames(GPL22945)))
pData(norm_GPL22945) <- pData(GPL22945)[sampleNames(norm_GPL22945),]


###################################################### GPL24120 ###################
# The same, but for the other platform

GPL24120 <- getGEO(filename = "./CKD_Data/GSE104948-GPL24120_series_matrix.txt.gz")
GEOACC_GPL24120 <- basename(as.character(pData(GPL24120)$supplementary_file))

# Read raw data
rawData_GPL24120 <- oligo::read.celfiles(fls[GEOACC_GPL24120])
sampleNames(rawData_GPL24120) <- sapply(sampleNames(rawData_GPL24120),function(sampleName) strsplit(sampleName,split="_")[[1]][1])

# RMA normalization
norm_GPL24120 <- oligo::rma(rawData_GPL24120,background=TRUE, normalize=TRUE)

# Update sample metadata
stopifnot(all(sampleNames(norm_GPL24120) %in% sampleNames(GPL24120)))
pData(norm_GPL24120) <- pData(GPL24120)[sampleNames(norm_GPL24120),]


# So both, norm_GPL22945 and norm_GPL24120 are two expression-sets with RMA-qt-normalized profiles
# splitted by platform from GSE104948 data set.

saveRDS(norm_GPL22945, file="./data/GSE104948_GPL570.rds")
saveRDS(norm_GPL24120, file="./data/GSE104948_GPL96.rds")
Collapse