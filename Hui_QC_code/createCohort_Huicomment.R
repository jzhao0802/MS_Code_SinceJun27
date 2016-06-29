rm(list=ls())
########################################################
#Hui: import plyr to use ldply function
########################################################
library(plyr)
library(dplyr)
library(caret)
# 001 create cohort table
# Comments by Hui on 6/29/2016

# 
main.wkDir <- "./"
#########################################################
#Hui: wk_dir does not exist, should be main.wkDir
########################################################
setwd(main.wkDir)

# main.inDir <- paste0(main.wkDir, '../01_Data/')
#Hui: data path
main.inDir <- "F:\\Jie\\MS\\01_Data\\" 

#Hui: time stamp
main.timeStamp <- as.character(Sys.time())
main.timeStamp <- gsub(":", ".", main.timeStamp)  # replace ":" by "."

###########################################################
#Hui: output path, shouldn't name it as Result because it is not the results, 
#     they are still our pre-model data. I suggest to create a subfolder in 01_data, 
#     and named as "descriptive_data"
##################################################################
main.outDir <- paste("F:/Hui/Project_2016/MS/01_Data/descriptive_data/",  main.timeStamp, "/", sep = '') 
dir.create(main.outDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")

#Hui: QC purpose
main.bCati <- F
main.bTest <- F

#Hui: data file
main.inFileNm <- "MS_decsupp_analset_20160614"
main.inFileExt <- ".csv"

######################################################
#Hui: cohort number, should be 5
######################################################
main.cohortLst <- 1:5 

#Hui: outcome names
main.outcomeLst <- c("edssprog"
                     , 'edssconf3'
                     , 'relapse_fu_any_01'
                     , 'relapse_or_prog'
                     , 'relapse_and_prog'
                     , 'relapse_or_conf')

#Hui: missing values
main.na_represents <- c('', 'NA', 'unknown', 'ambiguous')

#Hui: pre-difined catgorical variables
main.varDefCati <- c("idx_rx", 'gender', 'birth_region', 'init_symptom')

#Hui: threshold for merging
main.threshold4merge <- 0.1

#outDir <- main.outDir
#inDir <- main.inDir
#cohortLst <- main.cohortLst
#outcomeLst <- main.outcomeLst
#bCati <- main.bCati
#bTest <- main.bTest
#inFileNm <- main.inFileNm
#inFileExt <- main.inFileExt
#na_represents <- main.na_represents
#varDefCati <- main.varDefCati
#threshold4merge <- main.threshold4merge

##########################################################################################
#Hui: the whole process of createCohortTb is:
# 1. reading in original dta
# 2. split the original data into different cohorts, run loops on cohorts
# 3. remove some variables as discussed
# 4. create new variable pre_dmt for different years
# 5. convert variables:
#   5.1 identify variable types for each column
#   5.2 cut the continuous variables into buckets
#   5.3 for catigorical vari:
#       5.3.1 merge small levels in catigorical variable into 1 level with order
#       5.3.2 merge small levels in catigorical variable into 1 level without order
#       5.3.3 determine the type of variables, 
#           if transfer then all variables are catigorical(continous is cut into buckets)
#           if no-transfer then there are 2 type of variables: continuous and catigorical
#       5.3.4 1 variable with values of (0, -1, NA) need to be taken care of
#   5.4 convert NA to 999
#   5.5 convert catigorical and numeric variables to dummy variables
# 6. ouput cohort data by outcome
###########################################################################################


#Hui: create cohort data, prepare for descriptive table
#Hui: in total 5 cohorts, 6 outcomes, and 2 types = 60 tables
createCohortTb <- function(inDir, inFileNm, inFileExt, outDir
                           , cohortLst, outcomeLst, bCati, na_represents
                           , varDefCati, threshold4merge
                           #Hui: add bTest
                           , bTest){
  #########################################################
  #Hui: 1. reading data
  ##############################################
  dt <- read.table(paste0(inDir, inFileNm, inFileExt)
                   , sep=','
                   , header = T
                   , stringsAsFactors = F
                   , na.strings = na_represents)
  #Hui: QC 
  if(bTest == T){
    dt = dt[1:1000, ]
  }
  
  #Hui: convert variable names to lower case
  names(dt) <- tolower(names(dt))
  dim(dt) #[1] 6501  411
  
  # for a certain cohort, for those duplicated ptid , randomly select one line
  for(cohort in cohortLst){
    
    #######################################################################
    #Hui: 2. split the data into each cohort and run for loop on each cohort
    ######################################################################
    if(cohort == 5){
      dtCoh <- dt
    }else if(cohort %in% 1:4){
      dtCoh <- dt %>% filter(tblcoh==cohort)
    }else{
      stop("wrong input cohort index!\n")
    }
    
    ##################################################################################################
    #Hui: 3. remove some variables based on discussion between Jie and Lichao, see tracking file line ***
    ##################################################################################################
    dtCoh <- dtCoh %>%
      select(-idx_dt) %>%
      select(-firstdt) %>%
      group_by(new_pat_id) %>%
      do(sample_n(., 1)) %>%
#       select(-new_pat_id) %>%
      select(-tblcoh)
    
    #Hui: model variable name list
    varLst <- names(dtCoh)
    
    #Hui: define cohort number, using for output file name
    cohortNm <- ifelse(cohort==1, "BConti"
                       , ifelse(cohort==2, "B2B"
                                , ifelse(cohort== "B2Fir"
                                         , ifelse(cohort==4, "B2Sec"
                                                  , ifelse(cohort==5, "Cmp", step("wrong cohort index!\n"))))))
    
    ###################################################################################
    #Hui: 4. create new variable-pre_index_dmts, based on Lichao's email on June 24
    ###################################################################################
    # calculte pre index dmts in defferent yeas
    #Hui: 4.1 grep the name from varLst
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
    #Hui: 4.2
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
    
    #Hui 4.3 create pre_dmts for each year
    dtCoh$pre_dmts_1 <- apply(dtCoh_forDmts[, var_preDmt_1], 1, sum, na.rm=T)
    dtCoh$pre_dmts_2 <- apply(dtCoh_forDmts[, var_preDmt_2], 1, sum, na.rm=T)
    dtCoh$pre_dmts_3 <- apply(dtCoh_forDmts[, var_preDmt_3], 1, sum, na.rm=T)
    dtCoh$pre_dmts_4 <- apply(dtCoh_forDmts[, var_preDmt_4], 1, sum, na.rm=T)
    
    #Hui: QC
    dim(dtCoh) #[1] 383 413
    
    #Hui: final variable name list
    varLst_f1 <- names(dtCoh)
    
    ###############################################################
    #Hui: flag to determine whether there will be a transfermation
    #     this is actually working for output data naming, will be overwritten by line 203
    #     I suggest to use ifelse fucntion to determine the flag instead of overwritting
    ##########################################################
    
    ##########################################################
    #Hui: didn't define bTransf before
    ##########################################################
    flag <- ifelse(bTransf==T, 'withTransf', 'withoutTransf')
    #flag <- "withoutTransf"
    
    ##############################################################################
    #Hui: 5. convert the variables according to the type 
    #######################################################################
    if(bTransf==T){
      #flag <- "withTransf"
      #############################################################
      #Hui: 5.1 identify variable types for each column
      #########################################################
      varClfList <- varClassify(dtCoh)
      
      #####################################################
      #Hui: duplicated with line 59, varDefCati is an input of createCohort function
      ####################################################
      #varDefCati <- c("idx_rx", 'gender', 'birth_region', 'init_symptom')
      
      #Hui: variables to be merged
      var2merge <- setdiff(varDefCati, c(varClfList$naVars, varClfList$cansVars, varClfList$biVars))
      
      # new_pat_id should not be transformed
      #Hui: variables to be cut
      var2quartile <- setdiff(c(varClfList$catVars, varClfList$contVars), 
                              c(var2merge, varClfList$biVars, 'new_pat_id'))
      
      #Hui: variables to be convert to dummy
      var2quartileBnumeric <- setdiff(var2quartile, varClfList$charVars)
      
      #     temp <- with(dtCoh[, var2quartileBnumeric], cut(var2quartileBnumeric, 
      #                                     breaks=quantile(var2quartileBnumeric, probs=seq(0,1, by=0.25), na.rm=TRUE), 
      #                                     include.lowest=TRUE))
      
      #Hui: define data as a dataframe
      dtCoh <- as.data.frame(dtCoh)
      
      ######################################################
      #Hui: 5.2 cut the continuous variables into buckets
      ###########################################################
      dt2quartile <- as.data.frame(t(ldply(lapply(var2quartileBnumeric, function(var){
        
        varVct <- dtCoh[, var]
        rowQuartile <- as.character(cut(varVct
                                        , breaks=unique(quantile(varVct, probs=seq(0, 1, by=1/4), na.rm=T))
                                        , include.lowest = T))
        return(rowQuartile)
      }), quickdf)))
      
      names(dt2quartile) <- var2quartileBnumeric
      
      ##################################################################
      #Hui: 5.3 for catigorical variable:merge small levels into 1 level
      #################################################################
      var2quartileBchar <- setdiff(var2quartile, var2quartileBnumeric)
      
      #############################################################################
      #Hui: change the threshold to threshold4merge
      #  5.3.1 merge small levels in catigorical variable into 1 level with order
      #############################################################################
      dt2mergeGrad <- as.data.frame(t(ldply(lapply(var2quartileBchar
                                                   , function(var)merge4withGradCatiVars(var, dtCoh, threshold4merge))
                                            , quickdf)))
      names(dt2mergeGrad) <- var2quartileBchar
      
      ######################################################################
      #Hui: remove threshold becaste the argument in function is threshold4merge
      # 5.3.2 merge small levels in catigorical variable into 1 level without order
      #########################################################################
      #threshold <- 0.1
      dt2merge <- as.data.frame(t(ldply(lapply(var2merge
                                               , function(var)merge4CatiVars(var, dtCoh, threshold4merge ))
                                        , quickdf)))
      names(dt2merge) <- var2merge
      
      #Hui: transformed data
      dtCoh <- as.data.frame(
        cbind(dt2quartile
              , dt2mergeGrad
              , dt2merge
              , dtCoh[, setdiff(varLst_f1, c(var2quartileBnumeric, var2quartileBchar, var2merge))]))
      
    }
    
    # transfor all the charact variables into dummy using model.matrix
    #########################################################################
    #Hui: 5.3.3 determine the type of variables, 
    #           if transfer then all variables are catigorical(continous is cut into buckets)
    #           if no-transfer then there are 2 type of variables: continuous and catigorical
    #########################################################################
    varTypeLst <- getVarType(dt=dtCoh, varLst = colnames(dtCoh))
    charVars <- varTypeLst$charVars
    numVars <- varTypeLst$numVars

    # other numeric columns should be transformed into dummy
    ###########################################################
    #Hui: 5.3.4 1 variable with values of (0, -1, NA) need to be taken care of
    ###############################################################
    b2dummy <- sapply(dtCoh[,numVars], function(x){
      lvs <- unique(x)
      length(setdiff(lvs, c(0, 1, NA))) > 0
    })
    varNumB2dummy <- setdiff(numVars[b2dummy], "new_pat_id")
    
    charVars <- c(charVars, varNumB2dummy)
    
    dtCohChar <- dtCoh[, charVars]
    
    # dtCohChar2Fct <- sapply(as.data.frame(dtCoh[, charVars]), factor)
    
    #########################################################################
    #Hui: 5.4 convert NA to 999
    #########################################################################
    # before turn to dummy, replace NA using 999
    dtCohCharRepNA <- as.data.frame(t(ldply(lapply(charVars, function(var){
      vct <- dtCohChar[, var]
      char <- as.character(vct)
      char[is.na(char)] <- 999
      # fct <- as.factor(char)
      return(char)
    }), quickdf)))
    names(dtCohCharRepNA) <- charVars
    
    #Hui: convert to factor
    # turnto factor type
    dtCohChar2Fct <- as.data.frame(unclass(dtCohCharRepNA))
    
    #############################################################################
    #Hui: 5.5 convert catigorical and numeric variables to dummy variables
    ##############################################################################
    dtCohChar2Fct2Dummy <- getDummy(dtCohChar2Fct)
    dtCohFinal1 <- bind_cols(dtCohChar2Fct2Dummy
                            , dtCoh[, setdiff(varLst_f1, charVars)]) %>%
      as.data.frame(.)
    
    #############################################################################
    #Hui: 6. ouput cohort data by outcome
    #############################################################################
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








