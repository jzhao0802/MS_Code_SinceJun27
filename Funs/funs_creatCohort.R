
varClassify <- function(dt){
  levels_cnt <- sapply(dt, function(x)length(na.omit(unique(x))))
  varLst <- names(dt)
  naVars <- varLst[levels_cnt ==0]#2
  cansVars <- varLst[levels_cnt==1] #41
  biVars <- varLst[levels_cnt==2]   #277
  #   catVars <- varLst[levels <= 10 & levels >2]
  #   contVars <- varLst[levels > 10]
  catVars <- varLst[levels_cnt <= 20 & levels_cnt > 2]#84
  contVars <- varLst[levels_cnt > 20] #9
  # QC
  if(length(naVars) + length(cansVars) + length(biVars)  + length(catVars) + length(contVars) != dim(dt)[2]){
    stop("wrong!\n")
  }
  
  contVarsLvsLst <- lapply(contVars, function(var)table(dt[, var]))
  names(contVarsLvsLst) <- contVars
  contVarslvsNum <- lapply(contVarsLvsLst, length)
  catVarsLvsLst <- lapply(catVars, function(var)table(dt[, var]))
  names(catVarsLvsLst) <- catVars
  catVarslvsNum <- lapply(catVarsLvsLst, length)
  
  bNum <- sapply(dt, function(x)is.numeric(x))
  numVars <- varLst[bNum]
  charVars <- varLst[!bNum]
  varClf <- list(naVars=naVars, cansVars=cansVars, biVars=biVars, catVars=catVars, contVars=contVars, numVars=numVars, charVars=charVars)
  return(varClf)
}

merge4CatiVars <- function(var, dtCoh, threshold){
  vct <- dtCoh[, var]
  tb <- table(vct)
  prob <- prop.table(tb)
  totCnt <- sum(tb)
  refTb <- data.frame(cnt=tb, prob=prob)[, c(2,4)]
  names(refTb) <- c("cnt", 'prob')
  rownames(refTb) <- names(tb)
  
  refTb <- refTb[order(refTb$cnt), ]
  
  for(i in 1:nrow(refTb)){
    if(refTb$prob[1] < threshold){
      lvsB4merge <- rownames(refTb)[1:(1+1)]
      mergeLvs <-paste0(lvsB4merge, collapse = ' OR ')
      merge <- apply(refTb[1:(1+1),], 2, sum)
      refTb <- rbind(merge, refTb[-(1:(1+1)), ])
      rownames(refTb)[1] <- mergeLvs
      refTb <- refTb[order(refTb$cnt), ]
      vct <- ifelse(vct %in% lvsB4merge, mergeLvs, vct)
    }else{
      break
    }
  }
  return(vct)
}



merge4withGradCatiVars <- function(var, dtCoh, threshold){
  vct <- dtCoh[, var]
  tb <- table(vct)
  prob <- prop.table(tb)
  totCnt <- sum(tb)
  refTb <- data.frame(cnt=tb, prob=prob)[, c(2,4)]
  names(refTb) <- c("cnt", 'prob')
  rownames(refTb) <- names(tb)
  if(grepl('cranial|spinal', var, ignore.case = T)){
    refTb <- refTb[c(2:nrow(refTb), 1), ]
  }

  for(i in 1:nrow(refTb)){
    minIdx <- which(refTb$prob == min(refTb$prob))
    if(length(minIdx)>1){
      minIdx <- sample(minIdx, 1)
    }
    minProb <- refTb$prob[minIdx]
    if(minProb < threshold ){
      if(minIdx == 1){
        idx2merge <- c(1,2)
        lvsB4merge <- rownames(refTb)[idx2merge]
        mergeLvs <-paste0(lvsB4merge, collapse = ' OR ')
        merge <- apply(refTb[idx2merge,], 2, sum)
        refTb <- rbind(merge, refTb[-idx2merge, ])
        rownames(refTb)[1] <- mergeLvs
        
      }else if(minIdx == nrow(refTb)){
        idx2merge <- c((nrow(refTb)-1):nrow(refTb))
        lvsB4merge <- rownames(refTb)[idx2merge]
        mergeLvs <-paste0(lvsB4merge, collapse = ' OR ')
        merge <- apply(refTb[idx2merge,], 2, sum)
        refTb <- rbind(refTb[-idx2merge, ], merge)
        rownames(refTb)[nrow(refTb)] <- mergeLvs
        
      }else{
        leftRight <- refTb$prob[c(minIdx-1, minIdx+1)]
        idx2mergeAno <- c(minIdx-1, minIdx+1)[leftRight==min(leftRight)]
        idx2merge <- c(minIdx, idx2mergeAno)
        lvsB4merge <- rownames(refTb)[idx2merge]
        mergeLvs <-paste0(lvsB4merge, collapse = ' OR ')
        merge <- apply(refTb[idx2merge,], 2, sum)
        if(min(idx2merge)==1){
          refTb <- rbind(merge, refTb[-idx2merge, ])
          rownames(refTb)[1] <- mergeLvs

        }else{
          refTb <- rbind(refTb[min(idx2merge)-1,]
                         , merge
                         , refTb[-(1:max(idx2merge)),])
          rownames(refTb)[min(idx2merge)] <- mergeLvs
        }

      }
      
      vct <- ifelse(vct %in% lvsB4merge, mergeLvs, vct)
    }else{
      break
    }
  }
  return(vct)
}


getVarType <- function(dt, varLst){
  varClass <- sapply(dt, function(x)class(x))
  lgVars <- varLst[varClass=="logical"] #2 vars are all NA
  charVars <- varLst[varClass %in% c("character", "factor")] #5
  numVars <- varLst[varClass %in% c('numeric', 'integer')] #406
  typeList <- list(lgVars=lgVars, charVars=charVars, numVars=numVars)
  return(typeList)
}


getDummy <- function(temp_fct){
  options(na.action="na.pass")
  lvsCnt <- sapply(temp_fct, function(x){length(levels(x))})
  var1lvs <- names(lvsCnt[lvsCnt<2])
  var2dummy <- setdiff(names(temp_fct), var1lvs)
  dummy <- 
    model.matrix( ~ .
                  , data=temp_fct[, var2dummy]
                  , contrasts.arg = 
                    lapply(temp_fct[, var2dummy]
                           , contrasts, contrasts=FALSE)
                  , na.action=na.pass
    )[, -1]
  if(length(var1lvs)>0){
    feakRows <- 2
    feakDt <- temp_fct[1:feakRows, var1lvs]
    feakDt[,] <- 999
    feakDt2dummy <- rbind(feakDt, temp_fct[, var1lvs])
    feakDt2dummy <- as.data.frame(unclass(feakDt2dummy))
    dummy1lvs <- 
      model.matrix( ~ .
                    , data=feakDt2dummy
                    , contrasts.arg = 
                      lapply(feakDt2dummy
                             , contrasts, contrasts=FALSE)
                    , na.action=na.pass
      )[, -1]
    dummy1lvsRm999 <- dummy1lvs[-(1:feakRows), !grepl("999$", colnames(dummy1lvs))]
    dummyAll <- as.data.frame(cbind(dummy, dummy1lvsRm999))
    
  }else{
    dummyAll <- as.data.frame(dummy)
  }
  
  return(dummyAll)
}

createCohortTb <- function(inDir, inFileNm, inFileExt, outDir
                           , cohortLst, outcomeLst, bTransf, na_represents
                           , varDefCati, threshold, bTest, bQcMode){
  
  dt <- read.table(paste0(inDir, inFileNm, inFileExt)
                   , sep=','
                   , header = T
                   , stringsAsFactors = F
                   , na.strings = na_represents)
  cat("data readin successfully!\n")
  if(bTest == T){
    dt = dt[1:1000, ]
  }
  
  names(dt) <- tolower(names(dt))
  dim(dt) #[1] 6501  411
  
  # for a certain cohort, for those duplicated ptid , randomly select one line
  for(cohort in cohortLst){
    cat('cohort:', cohort, ' start!\n')
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
      select(-idxyr) %>%
      group_by(new_pat_id) %>%
      do(sample_n(., 1)) %>%
      #       select(-new_pat_id) %>%
      select(-tblcoh)
    varLst <- names(dtCoh)
    cat('\nline sample for duplicated patid!\n')
    
    if(bQcMode==T){
      if(nrow(dtCoh)!=length(unique((dtCoh$new_pat_id)))){
        stop("unique patient id select wrong!\n\n")
      }
    }
    cohortNm <- ifelse(cohort==1, "BConti"
                       , ifelse(cohort==2, "B2B"
                                , ifelse(cohort==3, "B2Fir"
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
    if(bQcMode == T){
      if(any(dtCoh$pre_dmts_1 < apply(dtCoh[, var_preDmt_1], 1, sum, na.rm=T))
         & all(dtCoh$pre_dmts_2 < apply(dtCoh[, var_preDmt_2], 1, sum, na.rm=T))
         & all(dtCoh$pre_dmts_3 < apply(dtCoh[, var_preDmt_3], 1, sum, na.rm=T))
         & all(dtCoh$pre_dmts_4 < apply(dtCoh[, var_preDmt_4], 1, sum, na.rm=T))
         ){
        stop("the pre Dmts number is wrong!\n\n")
      }
      rxLst <- c("rx_fing_", "rx_ga_"
                 , "rx_nat_", "rx_ext_"
                 , "rx_bet_", "rx_avo_"
                 , "rx_reb_", "rx_tecf_"
                 , "rx_teri_", "rx_alem_")
      for(yr in 1:4){
        qcFlag <- unlist(lapply(1:10, function(rxIdx){
          dtCohRx <- dtCoh %>% filter(idx_rx==rxIdx & paste0(rxLst, yr)[rxIdx]==yr)
          eval(
            parse(
              text=
                paste0("flag <- apply(dtCohRx[, var_preDmt_"
                       , yr
                       , "], 1, sum, na.rm=T) -1 == dtCohRx$pre_dmts_"
                       , yr)
              )
            )
          return(flag)
        }))
        if(any(!qcFlag)){
          stop("the preDmts number is wrong!\n\n")
        }
      }
      
    }
    cat('pre_dmts get successfully!\n')
    dim(dtCoh) #[1] 383 412
    varLst_f1 <- names(dtCoh)
    flag <- "withoutTransf"
    
    if(bQcMode==T){
      if(any(c("idx_dt", "firstdt", 'idxyr', 'tblcoh') %in% varLst_f1)){
        stop("the 4 variables, idx_dt, firstdt, idxyr, tblcoh, have not been removed completely!\n")
      }  
    }
    
    if(bTransf==T){
      cat('bTransf:', bTransf, '\n')
      
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
      cat('\nfor var2quartileBnumeric, brak into quartile successfully!\n')
      
      var2quartileBchar <- setdiff(var2quartile, var2quartileBnumeric)
      
      dt2mergeGrad <- as.data.frame(t(ldply(lapply(var2quartileBchar
                                                   , function(var)merge4withGradCatiVars(var, dtCoh, threshold))
                                            , quickdf)))
      names(dt2mergeGrad) <- var2quartileBchar
      cat("\nfor var2quartileBchar, mergeGrad successfully!\n")
      dt2merge <- as.data.frame(t(ldply(lapply(var2merge
                                               , function(var)merge4CatiVars(var, dtCoh, threshold ))
                                        , quickdf)))
      names(dt2merge) <- var2merge
      cat("\n for var2merge, mergeWithoutGrad successfully!\n")
      dtCoh <- as.data.frame(
        cbind(dt2quartile
              , dt2mergeGrad
              , dt2merge
              , dtCoh[, setdiff(varLst_f1, c(var2quartileBnumeric, var2quartileBchar, var2merge))]))
      
    }
    # transfor all the charact variables into dummy using model.matrix
    
    varTypeLst <- getVarType(dt=dtCoh, varLst = colnames(dtCoh))
    charVars <- varTypeLst$charVars
    if(bQcMode==T){
      type <- sapply(dtCoh[, charVars], function(vct)class(vct))
      if(any(!type %in% c("character", "factor"))){
        stop("the character and factor variable list is not right!\n\n")
      }
    }
    
    if(bTransf==T){
      numVars <- varTypeLst$numVars
      b2dummy <- sapply(dtCoh[,numVars], function(x){
        lvs <- unique(x)
        length(setdiff(lvs, c(0, 1, NA))) > 0 & sum(!is.na(lvs)) < 3
      })
      varNumB2dummy <- setdiff(numVars[b2dummy], "new_pat_id")
      
      charVars <- c(charVars, varNumB2dummy)
      cat('\n for bTransf==T, add varNumB2dummy into charVars wich will be dummy later!\n')
    }
    # other numeric columns should be transformed into dummy
    
    dtCohChar <- as.data.frame(dtCoh[, charVars])
    
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
    
    if(bQcMode==T){
      naCnt <- apply(apply(dtCohCharRepNA, 2, is.na), 2, sum)
      if(all(naCnt) != 0){
        stop("NAs have not been competely replaced in character columns!\n")
      }
    }
    # turnto factor type
    dtCohChar2Fct <- as.data.frame(unclass(dtCohCharRepNA))
    
    dtCohChar2Fct2Dummy <- getDummy(dtCohChar2Fct)
    
    if(bQcMode==T){
      otherLevExists <- sapply(dtCohChar2Fct2Dummy, function(vct)length(setdiff(unique(vct), c(0, 1)))>1)
      if(sum(otherLevExists)>1){
        stop("not all character variables are turned to 0 1 varaibles!\n")
      }
    }
    dtCohFinal1 <- bind_cols(dtCohChar2Fct2Dummy
                             , dtCoh[, setdiff(varLst_f1, charVars)]) %>%
      as.data.frame(.)
    cat("\ndummy final!\n")
    
#     re <- lapply(outcomeLst, function(outcome){
#       dtCohFinal1$response <- dtCohFinal1[, outcome]
#       # remove outcome varibles list
#       dtCohFinal <- dtCohFinal1[, -match(outcomeLst, names(dtCohFinal1))]
#       #       dtCoh$tblcoh <- NULL
#       write.table(dtCohFinal
#                   , paste0(outDir, 'dt_', cohortNm, '_', outcome, "_", flag, '.csv')
#                   , sep=','
#                   , row.names = F)
#       return("export cohort successfully!\n")
#     })
#     cat("\nexport final 6 tables successfully!\n")
    
    
    write.table(dtCohFinal1
                , paste0(outDir, 'dt_', cohortNm, '_', flag, '.csv')
                , sep=','
                , row.names=F
                , na=""
                )
    cat("\nexport final 1 tables with 5 response successfully!\n")
  }
  
#   return(re)
}
