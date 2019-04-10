#Define rename file function 

my.file.rename <- function(from, to) {
  todir <- dirname(to)
  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
  file.rename(from = from,  to = to)
}

#Define MoveFiles function

MoveFiles <- function(exper, platform, mat, wd){
  if(file.exists(file.path(wd, 'Documents', 'CKD_Data', exper, paste0(exper, '_RAW/', platform)))){
    lapply(lapply(colnames(mat), function(x){
      list.files(path = file.path(wd, 'Documents', 'CKD_Data', exper, paste0(exper, '_RAW')), pattern = x)}),
      function(x) {
        my.file.rename(from = file.path(wd, 'Documents', 'CKD_Data', exper, paste0(exper, '_RAW/', x)),
                       to = file.path(wd, 'Documents', 'CKD_Data', exper, paste0(exper, '_RAW/', platform, '/', x)))
      })
  } else {
    dir.create(file.path(wd,'Documents', 'CKD_Data', exper,paste0(exper, '_RAW/'), platform))
    lapply(lapply(colnames(mat), function(x){
      list.files(path = file.path(wd, 'Documents', 'CKD_Data',exper,paste0(exper, '_RAW')), pattern = x)}),
      function(x) {
        my.file.rename(from = file.path(wd,'Documents', 'CKD_Data',exper,paste0(exper, '_RAW/',x)),
                       to = file.path(wd, 'Documents', 'CKD_Data', exper,paste0(exper, '_RAW/',platform,'/',x)))
      })
  }
}


#Download libraries

library(GEOquery)
library(affy)
library(tidyverse)
library(Biobase)

wd <- file.path(getwd(),"..")
exper <- 'GSE104954'

#Download CEL files

iter <- 0
filePaths<- NULL
while(is.null(filePaths) && iter <= 100){
  iter <- iter +1
  try(
    filePaths <- getGEOSuppFiles(exper, baseDir = paste0(wd,'/Documents/CKD_Data'))
  )
}

#Extract CEL files from TAR file

if(file.exists(file.path(wd, '/Documents/CKD_Data', 'exper', paste0(exper, '_RAW')))) {
  untar(rownames(filePaths)[1], exdir = file.path(wd, "/Documents/CKD_Data", exper, paste0(exper,'_RAW')))
} else {dir.create(file.path(wd, 'CKD_Data', exper,paste0(exper, '_RAW')))
  untar(rownames(filePaths)[1], exdir = file.path(wd, "/Documents/CKD_Data", exper, paste0(exper,'_RAW')))
}


#Separating the two data sets

gse <- getGEO(exper, GSEMatrix = TRUE)
gse1 <- gse[[1]]; gse2 <- gse[[2]]
mat1 <- as.data.frame(exprs(gse[[1]])); mat2 <- as.data.frame(exprs(gse[[2]]))


#Define Platforms

platform1 <- 'GPL570'
platform2 <- 'GPL96'


#Assigning each platform (GPL96 & GPL570) with its respective data set

MoveFiles(exper, platform1, mat1, wd)       #GPL570
MoveFiles(exper, platform2, mat2, wd)       #GPL96
      

#Normalize data with RMA            
      
      cel_data1 <- ReadAffy(celfile.path = file.path(wd, 'Documents', 'CKD_Data', exper, paste0(exper, '_RAW'), platform1))
      cel_data2 <- ReadAffy(celfile.path = file.path(wd, 'Documents', 'CKD_Data', exper, paste0(exper, '_RAW'), platform2))
      
      eset1 <- rma(cel_data1)
      eset2 <- rma(cel_data2)
      
      
#create data frame for expression data
            
      mat1 <- as.data.frame(exprs(eset1))
      mat2 <- as.data.frame(exprs(eset2))

      
#create data frame for patient data
      
      gpl1 <- getGEO(platform1)
      GeneSymbol1 <- gpl1@dataTable@table$`Gene Symbol`
      gpl2 <- getGEO(platform2)
      GeneSymbol2 <- gpl2@dataTable@table$'Gene Symbol'
      
      
#To account for multiple probe measurements of the same gene, average probe measurements
      
      mat1 <- aggregate(mat1, by = list(GeneSymbol1), FUN = mean)
      colnames(mat1) <- gsub('_.+','',colnames(mat1))
      colnames(mat1)[1] <- 'GeneSymbol'
      
      mat2 <- aggregate(mat2, by = list(GeneSymbol2), FUN = mean)
      colnames(mat2) <- gsub('_.+','',colnames(mat2))
      colnames(mat2)[1] <- 'GeneSymbol'
      
