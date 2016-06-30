rm(list=ls())
library(plyr)
library(dplyr)
library(caret)
# 001 create cohort table

# 
main.wkDir <- "./"
setwd(main.wkDir)
source("./Funs/funs_creatCohort.R")
# main.inDir <- paste0(main.wkDir, '../01_Data/')
main.inDir <- "F:\\Jie\\MS\\01_Data\\"

main.timeStamp <- as.character(Sys.time())
main.timeStamp <- gsub(":", ".", main.timeStamp)  # replace ":" by "."
main.outDir <- paste("./03_Result/",  main.timeStamp, "/", sep = '')
dir.create(main.outDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
main.bTransf <- T
main.bTest <- F
main.inFileNm <- "MS_decsupp_analset_20160614"
main.inFileExt <- ".csv"

main.cohortLst <- 1:5

main.outcomeLst <- c("edssprog"
                     , 'edssconf3'
                     , 'relapse_fu_any_01'
                     , 'relapse_or_prog'
                     , 'relapse_and_prog'
                     , 'relapse_or_conf')

main.na_represents <- c('', 'NA', 'unknown', 'ambiguous')

main.varDefCati <- c("idx_rx", 'gender', 'birth_region', 'init_symptom')

main.threshold4merge <- 0.1

# outDir <- main.outDir
# inDir <- main.inDir
# cohortLst <- main.cohortLst
# outcomeLst <- main.outcomeLst
# bCati <- main.bCati
# bTest <- main.bTest
# inFileNm <- main.inFileNm
# inFileExt <- main.inFileExt
# na_represents <- main.na_represents
# varDefCati <- main.varDefCati
# threshold <- main.threshold4merge


temp_re <- createCohortTb(inDir=main.inDir
                          , inFileNm=main.inFileNm
                          , inFileExt=main.inFileExt
                          , outDir=main.outDir
                          , cohortLst=main.cohortLst
                          , outcomeLst=main.outcomeLst
                          , bTransf=F
                          , na_represents=main.na_represents
                          , varDefCati=main.varDefCati
                          , threshold=main.threshold4merge
                          , bTest=main.bTest
                          , bQcMode = T)





