
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
