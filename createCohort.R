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
main.outcomeLst <- c("edssprog"
                     , 'edssconf3'
                     , 'relapse_fuany_01'
                     , 'relapse_or_prog'
                     , 'relapse_and_prog'
                     , 'relapse_or_conf')
outDir <- main.outDir
inDir <- main.inDir
cohortLst <- main.cohortLst
outcomeLst <- main.outcomeLst
bCati <- main.bCati
bTest <- main.bTest
inFileNm <- main.inFileNm
inFileExt <- main.inFileExt

createCohortTb <- function(inDir, inFileNm, inFileExt, outDir, cohortLst, outcomeLst, bCati){
  na_represents <- c('', 'NA', 'unknown', 'ambiguous')
  
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
  
  # remove lines whose base_line_edss_score is missing and remove idxyr
  dt <- dt %>% 
    filter(!is.na(baseline_edss_score)) %>% 
    select(-idxyr)%>% 
    select(-idx_dt) %>% 
    select(-firstdt)
  varLst <- names(dt)
  # [1] 6501  408
  missingDataNum <- apply(
    apply(dt, 2, function(x){is.na(x)})
    , 
    2
    ,sum)
  
  levels_cnt <- sapply(dt, function(x)length(na.omit(unique(x))))
  
  naVars <- varLst[levels_cnt ==0]#2
  cansVars <- varLst[levels_cnt==1] #41
  biVars <- varLst[levels_cnt==2]   #277
#   catVars <- varLst[levels <= 10 & levels >2]
#   contVars <- varLst[levels > 10]
  catVars <- varLst[levels_cnt <= 20 & levels_cnt > 2]#84
  contVars <- varLst[levels_cnt > 20] #9
  # QC
  length(naVars) + length(cansVars) + length(biVars)  + length(catVars) + length(contVars) == dim(dt)[2]
  contVarsLvsLst <- lapply(contVars, function(var)table(dt[, var]))
  names(contVarsLvsLst) <- contVars
  contVarslvsNum <- lapply(contVarsLvsLst, length)
  catVarsLvsLst <- lapply(catVars, function(var)table(dt[, var]))
  names(catVarsLvsLst) <- catVars
  catVarslvsNum <- lapply(catVarsLvsLst, length)
  
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
      group_by(new_pat_id) %>%
      do(sample_n(., 1))
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
    
    
    varClfList <- varClassify(dtCoh)
    varDefCati <- c("idx_rx", 'gender', 'birth_region', 'init_symptom')
    var2merge <- setdiff(varDefCati, c(varClfList$naVars, varClfList$cansVars, varClfList$biVars))
    
    var2quartile <- setdiff(c(varClfList$catVars, varClfList$contVars), c(var2merge, varClfList$biVars))
    var2quartileBnumeric <- setdiff(var2quartile, varClfList$charVars)
    
#     temp <- with(dtCoh[, var2quartileBnumeric], cut(var2quartileBnumeric, 
#                                     breaks=quantile(var2quartileBnumeric, probs=seq(0,1, by=0.25), na.rm=TRUE), 
#                                     include.lowest=TRUE))
    dtCoh <- as.data.frame(dtCoh)
    
    dt2quartile <- t(ldply(lapply(var2quartileBnumeric, function(var){
      
      varVct <- dtCoh[, var]
      rowQuartile <- as.character(cut(varVct
                         , breaks=unique(quantile(varVct, probs=seq(0, 1, by=1/4), na.rm=T))
                         , include.lowest = T))
      return(rowQuartile)
    }), quickdf))
    
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
    
    for(outcome in outcomeLst){
      dtCoh$response <- dtCoh[, outcome]
      # remove outcome varibles list
      dtCoh <- dtCoh[, -outcomeLst]
      dtCoh$tblcoh <- NULL

    }
    
  }
  
}








