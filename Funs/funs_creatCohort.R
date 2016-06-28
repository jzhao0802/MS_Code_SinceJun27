
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
  return(dummyAll)
}

