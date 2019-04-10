#Download libraries

library(GEOquery)
library(affy)
library(tidyverse)
library(Biobase)

#Download CEL files

exp <- 'GSE104948'
iter <- 0
filePaths <- NULL
while(is.null(filePaths) && iter <= 100){
  iter <- iter + 1
  try(
    filePaths <- getGEOSuppFiles(exp, baseDir = paste0(wd,'/Data/CELFiles'))
    
#Extract CEL files from TAR file
    
    if(file.exists(file.path(wd, 'Data', 'CELFiles', exp, paste0(exp, '_RAW')))) {
      untar(rownames(filePaths)[1], exdir = file.path(wd, "Data", "CELFiles", exp, paste0(exp,'_RAW')))
    } else {
      dir.create(file.path(wd, 'Data', 'CELFiles', exp,paste0(exp, '_RAW')))
      untar(rownames(filePaths)[1], exdir = file.path(wd, "Data", "CELFiles", exp, paste0(exp,'_RAW')))
      
#Separating the two data sets
      
      iter <- 0
      gse <- NULL
      while(is.null(gse) && iter <= 100){
        iter <- iter + 1
        try(
          gse <- getGEO(exp, GSEMatrix = TRUE)
        )
      }
      gse1 <- gse[[1]]; gse2 <- gse[[2]]
      mat1 <- as.data.frame(exprs(gse[[1]])); mat2 <- as.data.frame(exprs(gse[[2]]))
      
      
#Assigning each platform (GPL96 & GPL570) with its respective data set
      
      
