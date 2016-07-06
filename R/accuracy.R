accuracy <- function(pred, obs, cf = 0.5){
  obs <- factor(obs,levels=0:1)
  if (!is.matrix(pred)) 
  	dim(pred) <- c(length(pred),1)
  npred <- ncol(pred)
  L <- matrix(0,6,npred)
  rownames(L) <- c('Sens', 'Spec','Accu', 'Prec', 'BalancedAccu', 'F1')
  colnames(L) <- paste0("pred",1:npred)
  if(length(cf) == 1) cf = rep(cf, npred)
  for (i in 1:npred)
  {
	  predi <- factor(as.integer(pred[,i] >= cf[i]),levels=0:1)
	  tab <- table(predi,obs)
	  TP = tab[2,2]
	  TN = tab[1,1]
	  FP = tab[2,1]
	  FN = tab[1,2]
	  Sens = ifelse((TP+FN) == 0, 0, TP/(TP+FN))
	  Spec = ifelse((FP+TN) == 0, 0, TN/(FP+TN))
	  Accu = (TP+TN) / (TP+TN+FP+FN)
	  Prec = ifelse((TP+FP) == 0, 0, TP/(TP+FP))
	  BalancedAccu = (Sens + Spec) / 2.0
	  F1 = ifelse((Sens+Prec) == 0, 0, 2 * (Sens*Prec)/(Sens+Prec))
  	  L[,i] <- c(Sens, Spec, Accu, Prec, BalancedAccu, F1)
  }
  return(L)
}


nested.accuracy <- function(pred, obs, indecies, cf = 0.5){
  obs <- factor(obs,levels=0:1)
  if (!is.matrix(pred)) 
    dim(pred) <- c(length(pred),1)
  if (!is.matrix(cf)) 
    dim(cf) <- c(length(cf),1)
  if(nrow(cf) == 1)
    cf = matrix(cf, length(indecies[[1]]), ncol(pred))
  npred <- ncol(pred)
  L <- matrix(0,6,npred)
  rownames(L) <- c('Sens', 'Spec','Accu', 'Prec', 'BalancedAccu', 'F1')
  colnames(L) <- paste0("pred",1:npred)
  if(length(cf) == 1) cf = rep(cf, npred)
  outer.folds = length(indecies[[1]])
  predi = matrix(NA, nrow = length(obs), ncol = npred)
  for (i in 1:npred)
  {
    for(j in 1:outer.folds){
      slice = outer.indecies[[i]][[j]][[1]]
      predi[slice,i] <- ifelse(pred[slice,i] >= cf[j,i], 1, 0)
    }
    tab <- table(factor(predi[,i],levels = 0:1), obs)
    TP = tab[2,2]
    TN = tab[1,1]
    FP = tab[2,1]
    FN = tab[1,2]
    Sens = ifelse((TP+FN) == 0, 0, TP/(TP+FN))
    Spec = ifelse((FP+TN) == 0, 0, TN/(FP+TN))
    Accu = (TP+TN) / (TP+TN+FP+FN)
    Prec = ifelse((TP+FP) == 0, 0, TP/(TP+FP))
    BalancedAccu = (Sens + Spec) / 2.0
    F1 = ifelse((Sens+Prec) == 0, 0, 2 * (Sens*Prec)/(Sens+Prec))
    L[,i] <- c(Sens, Spec, Accu, Prec, BalancedAccu, F1)
  }
  return(L)
}
