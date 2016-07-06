roc = function(probs, y, thresh.vec = seq(0.01, 0.99, by = 0.01)){
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
  

  opt.thresh.ba = thresh.vec[which.max(BalancedAccu)]
  opt.thresh.f1 = thresh.vec[which.max(F1)]
  distToCorner01 = sqrt((1 - tpr)^2 + (fpr)^2)
  opt.thresh.dist = thresh.vec[which.min(distToCorner01)]
  
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
reportResults <- function(pred.probs, y.labs, cutoff, method = 'equal prior'){
  cat(paste('prediction on test set accuracy with cutoff', paste(cutoff, collapse = ','), 'using method', method))
  cat('\n')
  L <- accuracy(pred.probs,y.labs, cutoff)
  print(L)
  cat('\n ============================================= \n')
  DF <- data.frame(cbind(apply(as.matrix(L), 1, mean), apply(as.matrix(L), 1, sd)), stringsAsFactors = F)
  rownames(DF) = rownames(L)
  colnames(DF) = c('mean', 'sd')
  print(DF)
  cat('\n ============================================= \n')
  
}


nested.reportResults <- function(pred.probs, y.labs, indecies, cutoff, method = 'equal prior'){
  cat(paste('prediction on test set accuracy with cutoff', paste(cutoff, collapse = ','), 'using method', method))
  cat('\n')

  L <- nested.accuracy(pred.probs,y.labs, indecies, cutoff)
  print(L)
  cat('\n ============================================= \n')
  DF <- data.frame(cbind(apply(as.matrix(L), 1, mean), apply(as.matrix(L), 1, sd)), stringsAsFactors = F)
  rownames(DF) = rownames(L)
  colnames(DF) = c('mean', 'sd')
  print(DF)
  cat('\n ============================================= \n')
  
}


