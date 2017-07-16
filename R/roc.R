roc = function(probs, y, thresh.vec = seq(0.01, 0.99, by = 0.1)){
  tpr = rep(0,length(thresh.vec))
  fpr = rep(0,length(thresh.vec))
  Sens = rep(0,length(thresh.vec))
  Spec = rep(0,length(thresh.vec))
  Accu = rep(0,length(thresh.vec))
  Prec = rep(0,length(thresh.vec))
  F1 = rep(0,length(thresh.vec))
  BalancedAccu = rep(0,length(thresh.vec))
  
  
  for(i in 1:length(thresh.vec)){
    obj <- accuracy(probs, y, cf = thresh.vec[i])
    Sens[i] <- obj[1]
    Spec[i] <- obj[2]
    Accu[i] <- obj[3]
    Prec[i] <- obj[4]
    BalancedAccu[i] <- obj[5]
    F1[i] <- obj[6]
    tpr[i] = Sens[i]
    fpr[i] = 1 - Spec[i]
  }
  

  ## there might be a range of threshold values that result to max accuracy.
  ## for example when probabilites are close to 1 and 0, every threshold value
  ## will give perfect predictions. In such cases, use the middle value.
  ##
  # opt.thresh.ba = thresh.vec[which.max(BalancedAccu)]
  # opt.thresh.f1 = thresh.vec[which.max(F1)]
  # distToCorner01 = sqrt((1 - tpr)^2 + (fpr)^2)
  # opt.thresh.dist = thresh.vec[which.min(distToCorner01)]
  # 
  distToCorner01 = sqrt((1 - tpr)^2 + (fpr)^2)
  m.ba <- max(BalancedAccu)
  m.f1 <- max(F1)
  m.d  <- min(distToCorner01)
  
  range.ba <- which(BalancedAccu == m.ba)
  range.f1 <- which(F1 == m.f1)
  range.d   <- which(distToCorner01 == m.d)
  
  opt.ba <- ceiling((range.ba[length(range.ba)] + range.ba[1]) / 2)
  opt.f1 <- ceiling((range.f1[length(range.f1)] + range.f1[1]) / 2)
  opt.d  <- ceiling((range.d[length(range.d)] + range.d[1]) / 2)
  
  opt.thresh.ba = thresh.vec[opt.ba]
  opt.thresh.f1 = thresh.vec[opt.f1]
  opt.thresh.dist = thresh.vec[opt.d]
  ##########
  ind = order(fpr)
  fpr <- fpr[ind]
  tpr <- tpr[ind]
  m <- length(fpr)

  Sens <- Sens[ind]
  Spec <- Spec[ind]
  Accu <- Accu[ind]
  Prec <- Prec[ind]
  BalancedAccu <- BalancedAccu[ind]
  F1 <- F1[ind]
  
  
  x = fpr
  y = tpr
  xp <- c(x, x[m:1])
  yp <- c(numeric(m), y[m:1])
  n <- 2 * m
  p1 <- sum(xp[1:(n - 1)] * yp[2:n]) + xp[n] * yp[1]
  p2 <- sum(xp[2:n] * yp[1:(n - 1)]) + xp[1] * yp[n]
  AUC <- (0.5 * (p1 - p2))

  r = list(TPR = tpr, FPR = fpr, AUC = AUC, opt.thresh.dist = opt.thresh.dist,
           opt.thresh.f1 = opt.thresh.f1, opt.thresh.ba = opt.thresh.ba)

  return(r)
}

## Accuracy reports
reportResults <- function(pred.probs, y.labs, cutoff, method = 'equal prior', verbose = TRUE){
  if(verbose){
    cat(paste('prediction on test set accuracy with cutoff', paste(cutoff, collapse = ','), 'using method', method))
    cat('\n')
  }
  L <- accuracy(pred.probs,y.labs, cutoff)
  if(verbose){
    print(L)
    cat('\n ============================================= \n')
  }
  DF <- data.frame(cbind(apply(as.matrix(L), 1, mean), apply(as.matrix(L), 1, sd)), stringsAsFactors = F)
  rownames(DF) = rownames(L)
  colnames(DF) = c('mean', 'sd')
  
  # Compute AUC and Roc Curve
  if(!is.matrix(pred.probs)){
    pred.probs = matrix(pred.probs, ncol = 1)
  }
  
  ROCs <- apply(pred.probs, 2, function(x) roc(x, y.labs))
  
  all.AUCs <- unlist(lapply(ROCs, function(x) x$AUC))
  AUC.mean <- mean(all.AUCs)
  AUC.sd   <- sd(all.AUCs)
  best.AUC <- which.max(all.AUCs)
  ROC <- data.frame(FPR = ROCs[[best.AUC]]$FPR, TPR = ROCs[[best.AUC]]$TPR, stringsAsFactors = F) 
  
  DF <- rbind(DF, c(AUC.mean, AUC.sd))
  rownames(DF)[nrow(DF)] = 'AUC'
  
  if(verbose){
    print(DF)
    cat('\n ============================================= \n')
  }
  
  return(list(DF = DF, ROC = ROC))
}


nested.reportResults <- function(pred.probs, y.labs, indecies, cutoff, method = 'equal prior', verbose = TRUE){
  if(verbose){
    cat(paste('prediction on test set accuracy with cutoff', paste(cutoff, collapse = ','), 'using method', method))
    cat('\n')
  }

  L <- nested.accuracy(pred.probs,y.labs, indecies, cutoff)
  if(verbose){
    print(L)
    cat('\n ============================================= \n')
  }
  DF <- data.frame(cbind(apply(as.matrix(L), 1, mean), apply(as.matrix(L), 1, sd)), stringsAsFactors = F)
  rownames(DF) = rownames(L)
  colnames(DF) = c('mean', 'sd')
  
  ROCs <- apply(pred.probs, 2, function(x) roc(x, y.labs))
  
  all.AUCs <- unlist(lapply(ROCs, function(x) x$AUC))
  AUC.mean <- mean(all.AUCs)
  AUC.sd   <- sd(all.AUCs)
  best.AUC <- which.max(all.AUCs)
  ROC <- data.frame(FPR = ROCs[[best.AUC]]$FPR, TPR = ROCs[[best.AUC]]$TPR, stringsAsFactors = F) 
  
  DF <- rbind(DF, c(AUC.mean, AUC.sd))
  rownames(DF)[nrow(DF)] = 'AUC'
  
  if(verbose){
    print(DF)
    cat('\n ============================================= \n')
  }
  
  return(list(DF = DF, ROC = ROC))
}

## Accuracy reports
reportNovelResults <- function(pred.probs, cutoff, method = 'equal prior', verbose = TRUE){
  if(verbose){
    cat(paste('prediction of test set labels with cutoff', paste(cutoff, collapse = ','), 'using method', method))
    cat('\n')
  }
  
  if(ncol(pred.probs) > 1){
    if(length(cutoff) == 1) cutoff = rep(cutoff, ncol(pred.probs))
    labs = matrix(-2, nrow = nrow(pred.probs), ncol = ncol(pred.probs))
    for(i in 1:ncol(labs)){
      labs[,i] = ifelse(pred.probs[,i] >= cutoff[i], 1, 0)
    }
    
    pred.probs <- apply(labs, 1, mean)
    labs <- ifelse(pred.probs >= 0.5, 1, 0)
  }else{
    labs <- ifelse(pred.probs >= cutoff, 1, 0)
  }
  
  
  
  return(labs)
}

