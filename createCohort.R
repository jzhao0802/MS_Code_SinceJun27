rm(list=ls())
library(dplyr)
library(caret)
# 001 create cohort table

# 
main.wkDir <- "./"
setwd(main.wk_dir)

# main.inDir <- paste0(main.wkDir, '../01_Data/')
main.inDir <- "F:\\Jie\\MS\\01_Data\\"

main.timeStamp <- as.character(Sys.time())
main.timeStamp <- gsub(":", ".", main.timeStamp)  # replace ":" by "."
main.outDir <- paste("./03_Result/",  main.timeStamp, "/", sep = '')
dir.create(main.outDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
main.bCati <- F
main.bTest <- F
main.inFileNm <- "MS_decsupp_analset_20160614"
main.inFileExt <- ".csv"

main.cohortLst <- 1:6

main.outcomeLst <- c("edssprog"
                     , 'edssconf3'
                     , 'relapse_fu_any_01'
                     , 'relapse_or_prog'
                     , 'relapse_and_prog'
                     , 'relapse_or_conf')

main.na_represents <- c('', 'NA', 'unknown', 'ambiguous')

main.varDefCati <- c("idx_rx", 'gender', 'birth_region', 'init_symptom')

main.threshold4merge <- 0.1

outDir <- main.outDir
inDir <- main.inDir
cohortLst <- main.cohortLst
outcomeLst <- main.outcomeLst
bCati <- main.bCati
bTest <- main.bTest
inFileNm <- main.inFileNm
inFileExt <- main.inFileExt
na_represents <- main.na_represents
varDefCati <- main.varDefCati
threshold4merge <- main.threshold4merge

createCohortTb <- function(inDir, inFileNm, inFileExt, outDir
                           , cohortLst, outcomeLst, bCati, na_represents
                           , varDefCati, threshold4merge){

  dt <- read.table(paste0(inDir, inFileNm, inFileExt)
                   , sep=','
                   , header = T
                   , stringsAsFactors = F
                   , na.strings = na_represents)
  if(bTest == T){
    dt = dt[1:1000, ]
  }
  
  names(dt) <- tolower(names(dt))
  dim(dt) #[1] 6501  411
  
  # for a certain cohort, for those duplicated ptid , randomly select one line
  for(cohort in cohortLst){
    if(cohort == 5){
      dtCoh <- dt
    }else if(cohort %in% 1:4){
      dtCoh <- dt %>% filter(tblcoh==cohort)
    }else{
      stop("wrong input cohort index!\n")
    }
    dtCoh <- dtCoh %>%
      select(-idx_dt) %>%
      select(-firstdt) %>%
      group_by(new_pat_id) %>%
      do(sample_n(., 1)) %>%
#       select(-new_pat_id) %>%
      select(-tblcoh)
    varLst <- names(dtCoh)
    
    cohortNm <- ifelse(cohort==1, "BConti"
                       , ifelse(cohort==2, "B2B"
                                , ifelse(cohort== "B2Fir"
                                         , ifelse(cohort==4, "B2Sec"
                                                  , ifelse(cohort==5, "Cmp", step("wrong cohort index!\n"))))))
    
    # calculte pre index dmts in defferent yeas
    var_preDmt_1 <- grep('^rx_(fing|ga|nat|ext|avo|reb|bet|tecf|teri|alem)_1'
                         , varLst
                         , ignore.case = T
                         , value = T)
    
    # calculte pre index dmts in defferent yeas
    var_preDmt_2 <- grep('^rx_(fing|ga|nat|ext|avo|reb|bet|tecf|teri|alem)_2'
                         , varLst
                         , ignore.case = T
                         , value = T)
    
    # calculte pre index dmts in defferent yeas
    var_preDmt_3 <- grep('^rx_(fing|ga|nat|ext|avo|reb|bet|tecf|teri|alem)_3'
                         , varLst
                         , ignore.case = T
                         , value = T)
    
    # calculte pre index dmts in defferent yeas
    var_preDmt_4 <- grep('^rx_(fing|ga|nat|ext|avo|reb|bet|tecf|teri|alem)_4'
                         , varLst
                         , ignore.case = T
                         , value = T)
    # make sure that the index DMT type is removed from this counting
    dtCoh_forDmts <- ldply(lapply(1:nrow(dtCoh), function(irow){
      row <- dtCoh[irow, ]
      if(row$idx_rx == 1){
        row[, grep("^rx_fing_\\d$", varLst, ignore.case = T, value=T)] = 0
      }else if(row$idx_rx == 2){
        row[, grep("^rx_ga_\\d$", varLst, ignore.case = T, value=T)] = 0
      }else if(row$idx_rx == 3){
        row[, grep("^rx_nat_\\d$", varLst, ignore.case = T, value=T)] = 0
      }else if(row$idx_rx == 4){
        row[, grep("^rx_ext_\\d$", varLst, ignore.case = T, value=T)] = 0
      }else if(row$idx_rx == 5){
        row[, grep("^rx_bet_\\d$", varLst, ignore.case = T, value=T)] = 0
      }else if(row$idx_rx == 6){
        row[, grep("^rx_avo_\\d$", varLst, ignore.case = T, value=T)] = 0
      }else if(row$idx_rx == 7){
        row[, grep("^rx_reb_\\d$", varLst, ignore.case = T, value=T)] = 0
      }else if(row$idx_rx == 10){
        row[, grep("^rx_tecf_\\d$", varLst, ignore.case = T, value=T)] = 0
      }else if(row$idx_rx == 11){
        row[, grep("^rx_teri_\\d$", varLst, ignore.case = T, value=T)] = 0
      }else if(row$idx_rx == 13){
        row[, grep("^rx_alem_\\d$", varLst, ignore.case = T, value=T)] = 0
      }
      return(row)
    }), quickdf)
    
    dtCoh$pre_dmts_1 <- apply(dtCoh_forDmts[, var_preDmt_1], 1, sum, na.rm=T)
    dtCoh$pre_dmts_2 <- apply(dtCoh_forDmts[, var_preDmt_2], 1, sum, na.rm=T)
    dtCoh$pre_dmts_3 <- apply(dtCoh_forDmts[, var_preDmt_3], 1, sum, na.rm=T)
    dtCoh$pre_dmts_4 <- apply(dtCoh_forDmts[, var_preDmt_4], 1, sum, na.rm=T)
    
    dim(dtCoh) #[1] 383 413
    varLst_f1 <- names(dtCoh)
    flag <- "withoutTransf"
    if(bTransf==T){
      flag <- "withTransf"
      varClfList <- varClassify(dtCoh)
      varDefCati <- c("idx_rx", 'gender', 'birth_region', 'init_symptom')
      var2merge <- setdiff(varDefCati, c(varClfList$naVars, varClfList$cansVars, varClfList$biVars))
      
      # new_pat_id should not be transformed
      var2quartile <- setdiff(c(varClfList$catVars, varClfList$contVars), c(var2merge, varClfList$biVars, 'new_pat_id'))
      var2quartileBnumeric <- setdiff(var2quartile, varClfList$charVars)
      
      #     temp <- with(dtCoh[, var2quartileBnumeric], cut(var2quartileBnumeric, 
      #                                     breaks=quantile(var2quartileBnumeric, probs=seq(0,1, by=0.25), na.rm=TRUE), 
      #                                     include.lowest=TRUE))
      dtCoh <- as.data.frame(dtCoh)
      
      dt2quartile <- as.data.frame(t(ldply(lapply(var2quartileBnumeric, function(var){
        
        varVct <- dtCoh[, var]
        rowQuartile <- as.character(cut(varVct
                                        , breaks=unique(quantile(varVct, probs=seq(0, 1, by=1/4), na.rm=T))
                                        , include.lowest = T))
        return(rowQuartile)
      }), quickdf)))
      
      names(dt2quartile) <- var2quartileBnumeric
      
      var2quartileBchar <- setdiff(var2quartile, var2quartileBnumeric)
      
      dt2mergeGrad <- as.data.frame(t(ldply(lapply(var2quartileBchar
                                                   , function(var)merge4withGradCatiVars(var, dtCoh, threshold))
                                            , quickdf)))
      names(dt2mergeGrad) <- var2quartileBchar
      
      threshold <- 0.1
      dt2merge <- as.data.frame(t(ldply(lapply(var2merge
                                               , function(var)merge4CatiVars(var, dtCoh, threshold ))
                                        , quickdf)))
      names(dt2merge) <- var2merge
      
      dtCoh <- as.data.frame(
        cbind(dt2quartile
              , dt2mergeGrad
              , dt2merge
              , dtCoh[, setdiff(varLst_f1, c(var2quartileBnumeric, var2quartileBchar, var2merge))]))
      
    }
    # transfor all the charact variables into dummy using model.matrix
    
    varTypeLst <- getVarType(dt=dtCoh, varLst = colnames(dtCoh))
    charVars <- varTypeLst$charVars
    numVars <- varTypeLst$numVars
    # other numeric columns should be transformed into dummy
    b2dummy <- sapply(dtCoh[,numVars], function(x){
      lvs <- unique(x)
      length(setdiff(lvs, c(0, 1, NA))) > 0
    })
    varNumB2dummy <- setdiff(numVars[b2dummy], "new_pat_id")
    
    charVars <- c(charVars, varNumB2dummy)
    
    dtCohChar <- dtCoh[, charVars]
    
    # dtCohChar2Fct <- sapply(as.data.frame(dtCoh[, charVars]), factor)
    
    # before turn to dummy, replace NA using 999
    dtCohCharRepNA <- as.data.frame(t(ldply(lapply(charVars, function(var){
      vct <- dtCohChar[, var]
      char <- as.character(vct)
      char[is.na(char)] <- 999
      # fct <- as.factor(char)
      return(char)
    }), quickdf)))
    names(dtCohCharRepNA) <- charVars
    # turnto factor type
    dtCohChar2Fct <- as.data.frame(unclass(dtCohCharRepNA))
    
    dtCohChar2Fct2Dummy <- getDummy(dtCohChar2Fct)
    dtCohFinal1 <- bind_cols(dtCohChar2Fct2Dummy
                            , dtCoh[, setdiff(varLst_f1, charVars)]) %>%
      as.data.frame(.)
    
    re <- lapply(outcomeLst, function(outcome){
      dtCohFinal1$response <- dtCohFinal1[, outcome]
      # remove outcome varibles list
      dtCohFinal <- dtCohFinal1[, -match(outcomeLst, names(dtCohFinal1))]
      #       dtCoh$tblcoh <- NULL
      write.table(dtCohFinal
                  , paste0(outDir, 'dt_', cohortNm, '_', outcome, "_", flag, '.csv')
                  , sep=','
                  , row.names = F)
      return("export cohort successfully!\n")
    })
    
    
  }
  
  return(re)
}








