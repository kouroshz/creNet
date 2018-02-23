generatePriors <- function(ents, rels, evidence, verbose=TRUE)
{
  
  # ents and rels are assumed to *only* contain genes for which expression levels are available 
  # Use e.g. matchDataToKb to guarantee it
  
  ## For each hypothesis, identify the children and non-children 
  ## and compute their predicted values from rels and observed values from evidence
  
  ## Predicted values for children
  src.uid <- factor(rels$srcuid, levels = unique(rels$srcuid))
  child.sgn <- numeric(nrow(rels)) # "conflict" --> 0 
  child.sgn[rels$type == "increase"] = 1
  child.sgn[rels$type == "decrease"] = -1
  child.sgn <- factor(child.sgn, levels = c(-1,0,1))
  child.sgn <- split(child.sgn,src.uid)
  
  ## Observed values for children 
  child.uid <- split(rels$trguid,src.uid)
  child.val <- lapply(child.uid, getGeneVals, evidence=evidence)
  
  ## Observed values for non-children
  non.child.uid <- lapply(child.uid, function(x) setdiff(ents$uid,x))
  non.child.val <- lapply(non.child.uid, getGeneVals, evidence=evidence)
  
  
  u.hyps <- as.integer(levels(src.uid)) # unique hypotheses
  nhyps <- length(u.hyps) # number of unique hypotheses
  
  if (verbose == TRUE)
    cat("\n Computing pvalues")
  results = data.frame(matrix(0, nrow  = 2 * length(u.hyps), ncol = 12), stringsAsFactors = F)
  colnames(results) = c('uid', 'name', 'regulation', 'correct.pred', 'incorrect.pred', 'score',
                        'total.reachable', 'significant.reachable', 'total.ambiguous', 'significant.ambiguous',
                        'unknow', 'ternary.pval')
  
  if (verbose == TRUE)
    cat("\n Total number of hypothesis to consider:", length(u.hyps))
  

  for(p.s in 1:length(u.hyps)){
    #cat('.')
    results[(2*(p.s-1)+1), 1] = u.hyps[p.s]
    results[(2*p.s), 1]       = u.hyps[p.s]
    results[(2*(p.s-1)+1), 2] = ents$name[which(ents$uid == u.hyps[p.s])]
    results[(2*p.s), 2]       = ents$name[which(ents$uid == u.hyps[p.s])]
    results[(2*(p.s-1)+1), 3] = 'up'
    results[(2*p.s), 3]       = 'down'
    
    npp = length(which(child.sgn[[p.s]] == 1 &  child.val[[p.s]] == 1))
    npm = length(which(child.sgn[[p.s]] == 1 &  child.val[[p.s]] == -1))
    npz = length(which(child.sgn[[p.s]] == 1 &  child.val[[p.s]] == 0))
    
    nmp = length(which(child.sgn[[p.s]] == -1 &  child.val[[p.s]] == 1))
    nmm = length(which(child.sgn[[p.s]] == -1 &  child.val[[p.s]] == -1))
    nmz = length(which(child.sgn[[p.s]] == -1 &  child.val[[p.s]] == 0))
    
    nrp = length(which(child.sgn[[p.s]] == 0 &  child.val[[p.s]] == 1))
    nrm = length(which(child.sgn[[p.s]] == 0 &  child.val[[p.s]] == -1))
    nrz = length(which(child.sgn[[p.s]] == 0 &  child.val[[p.s]] == 0))
    
    nzp = length(which(non.child.val[[p.s]] == 1))
    nzm = length(which(non.child.val[[p.s]] == -1))
    nzz = length(which(non.child.val[[p.s]] == 0))
    
    #pvalq = runCRE(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz, method = 'Quaternary')
    pvalt = runCRE(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz, method = 'Ternary')
    #pvale = runCRE(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz, method = 'Enrichment')
    ##pvalf = runCRE(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz, method = 'Fisher')
    
    
    qPlus  = npp + npm + npz
    qMinus = nmp + nmm + nmz
    qR     = nrp + nrm + nrz
    qZero  = nzp + nzm + nzz
    
    results[(2*(p.s-1)+1), 4]  = npp + nmm
    results[(2*(p.s-1)+1), 5]  = npm + nmp
    results[(2*(p.s-1)+1), 6]  = npp + nmm - (npm + nmp)
    results[(2*(p.s-1)+1), 7]  = qPlus + qMinus + qR
    results[(2*(p.s-1)+1), 8]  = npp + npm + nmp + nmm + nrp + nrm
    results[(2*(p.s-1)+1), 9]  = qR
    results[(2*(p.s-1)+1), 10] = nrp + nrm
    results[(2*(p.s-1)+1), 11] = qZero
    results[(2*(p.s-1)+1), 12] = pvalt$pval.up
    ##results[(2*(p.s-1)+1), 13] = pvalt$pval.up
    ##results[(2*(p.s-1)+1), 14] = pvale$pval.up
    
    results[(2*p.s), 4]  = nmp + npm
    results[(2*p.s), 5]  = npp + nmm
    results[(2*p.s), 6]  = nmp + npm - (npp + nmm)
    results[(2*p.s), 7]  = qPlus + qMinus + qR
    results[(2*p.s), 8]  = npp + npm + nmp + nmm + nrp + nrm
    results[(2*p.s), 9]  = qR
    results[(2*p.s), 10] = nrp + nrm
    results[(2*p.s), 11] = qZero
    results[(2*p.s), 12] = pvalt$pval.down
    ##results[(2*p.s), 13] = pvalt$pval.down
    ##results[(2*p.s), 14] = pvale$pval.down
  }
  
  
  ##results[,15] = fdrtool(results[,12], "pvalue", F, F, F, "fndr")$q
  ##results[,16] = fdrtool(results[,13], "pvalue", F, F, F, "fndr")$q
  ##results[,17] = fdrtool(results[,14], "pvalue", F, F, F, "fndr")$q
  
  #results = results[order(as.numeric(results$quaternary.pval)), ]
  #rownames(results) = 1:nrow(results)
  
  return(results)
  
}

generatePriorsOrig <- function(ents, rels, evidence, verbose=TRUE)
{
  
  ents.mRNA = ents[which(ents$type == 'mRNA'), ]
  
  ## For each hypothesis, identify the children and non-children and thier evidence values
  u.hyps = unique(rels$srcuid)
  child.uid = lapply(u.hyps, function(x) rels$trguid[which(rels$srcuid == x)])
  child.sgn = lapply(u.hyps, function(x) ifelse(rels$type[which(rels$srcuid == x)] == 'increase',
                                                1, ifelse(rels$type[which(rels$srcuid == x)] == 'decrease', -1, 0)))
  
  child.val = lapply(child.uid, function(x) getGeneVals(x, evidence))
  
  non.child.uid = lapply(child.uid, function(x) unique(ents.mRNA$uid[which(!(ents.mRNA$uid %in% x))]))
  non.child.val = lapply(non.child.uid, function(x) getGeneVals(x, evidence))
  
  ## Get the data slices corresponding to each hypothesis
  child.id = lapply(child.uid, function(x) as.numeric(ents.mRNA$id[match(x,ents.mRNA$uid)])) ## to get the id
  
  if (verbose == TRUE)
    cat("\n Computing pvalues")
  results = data.frame(matrix(0, nrow  = 2 * length(u.hyps), ncol = 12), stringsAsFactors = F)
  colnames(results) = c('uid', 'name', 'regulation', 'correct.pred', 'incorrect.pred', 'score',
                        'total.reachable', 'significant.reachable', 'total.ambiguous', 'significant.ambiguous',
                        'unknow', 'ternary.pval')
  
  if (verbose == TRUE)
    cat("\n Total number of hypothesis to consider:", length(u.hyps))
  
  
  for(p.s in 1:length(u.hyps)){
    #cat('.')
    results[(2*(p.s-1)+1), 1] = u.hyps[p.s]
    results[(2*p.s), 1]       = u.hyps[p.s]
    results[(2*(p.s-1)+1), 2] = ents$name[which(ents$uid == u.hyps[p.s])]
    results[(2*p.s), 2]       = ents$name[which(ents$uid == u.hyps[p.s])]
    results[(2*(p.s-1)+1), 3] = 'up'
    results[(2*p.s), 3]       = 'down'
    
    npp = length(which(child.sgn[[p.s]] == 1 &  child.val[[p.s]] == 1))
    npm = length(which(child.sgn[[p.s]] == 1 &  child.val[[p.s]] == -1))
    npz = length(which(child.sgn[[p.s]] == 1 &  child.val[[p.s]] == 0))
    
    nmp = length(which(child.sgn[[p.s]] == -1 &  child.val[[p.s]] == 1))
    nmm = length(which(child.sgn[[p.s]] == -1 &  child.val[[p.s]] == -1))
    nmz = length(which(child.sgn[[p.s]] == -1 &  child.val[[p.s]] == 0))
    
    nrp = length(which(child.sgn[[p.s]] == 0 &  child.val[[p.s]] == 1))
    nrm = length(which(child.sgn[[p.s]] == 0 &  child.val[[p.s]] == -1))
    nrz = length(which(child.sgn[[p.s]] == 0 &  child.val[[p.s]] == 0))
    
    nzp = length(which(non.child.val[[p.s]] == 1))
    nzm = length(which(non.child.val[[p.s]] == -1))
    nzz = length(which(non.child.val[[p.s]] == 0))
    
    #pvalq = runCRE(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz, method = 'Quaternary')
    pvalt = runCRE(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz, method = 'Ternary')
    #pvale = runCRE(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz, method = 'Enrichment')
    ##pvalf = runCRE(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz, method = 'Fisher')
    
    
    qPlus  = npp + npm + npz
    qMinus = nmp + nmm + nmz
    qR     = nrp + nrm + nrz
    qZero  = nzp + nzm + nzz
    
    results[(2*(p.s-1)+1), 4]  = npp + nmm
    results[(2*(p.s-1)+1), 5]  = npm + nmp
    results[(2*(p.s-1)+1), 6]  = npp + nmm - (npm + nmp)
    results[(2*(p.s-1)+1), 7]  = qPlus + qMinus + qR
    results[(2*(p.s-1)+1), 8]  = npp + npm + nmp + nmm + nrp + nrm
    results[(2*(p.s-1)+1), 9]  = qR
    results[(2*(p.s-1)+1), 10] = nrp + nrm
    results[(2*(p.s-1)+1), 11] = qZero
    results[(2*(p.s-1)+1), 12] = pvalt$pval.up
    ##results[(2*(p.s-1)+1), 13] = pvalt$pval.up
    ##results[(2*(p.s-1)+1), 14] = pvale$pval.up
    
    results[(2*p.s), 4]  = nmp + npm
    results[(2*p.s), 5]  = npp + nmm
    results[(2*p.s), 6]  = nmp + npm - (npp + nmm)
    results[(2*p.s), 7]  = qPlus + qMinus + qR
    results[(2*p.s), 8]  = npp + npm + nmp + nmm + nrp + nrm
    results[(2*p.s), 9]  = qR
    results[(2*p.s), 10] = nrp + nrm
    results[(2*p.s), 11] = qZero
    results[(2*p.s), 12] = pvalt$pval.down
    ##results[(2*p.s), 13] = pvalt$pval.down
    ##results[(2*p.s), 14] = pvale$pval.down
  }
  
  
  ##results[,15] = fdrtool(results[,12], "pvalue", F, F, F, "fndr")$q
  ##results[,16] = fdrtool(results[,13], "pvalue", F, F, F, "fndr")$q
  ##results[,17] = fdrtool(results[,14], "pvalue", F, F, F, "fndr")$q
  
  #results = results[order(as.numeric(results$quaternary.pval)), ]
  #rownames(results) = 1:nrow(results)
  
  return(results)
}

generatePriorsDavid <- function(ents, rels, evidence, verbose=TRUE)
{

	# ents and rels are assumed to *only* contain genes for which expression levels are available 
	# Use e.g. matchDataToKb to guarantee it
	
	## For each hypothesis, identify the children and non-children 
	## and compute their predicted values from rels and observed values from evidence
	
	## Predicted values for children
	src.uid <- factor(rels$srcuid, levels = unique(rels$srcuid))
	child.sgn <- numeric(nrow(rels)) # "conflict" --> 0 
	child.sgn[rels$type == "increase"] = 1
	child.sgn[rels$type == "decrease"] = -1
	child.sgn <- factor(child.sgn, levels = c(-1,0,1))
	child.sgn <- split(child.sgn,src.uid)

	## Observed values for children 
	child.uid <- split(rels$trguid,src.uid)
	child.val <- lapply(child.uid, getGeneVals, evidence=evidence)

	## Observed values for non-children
	non.child.uid <- lapply(child.uid, function(x) setdiff(ents$uid,x))
	non.child.val <- lapply(non.child.uid, getGeneVals, evidence=evidence)
	
			
	if (verbose)
		cat("\n Computing p-values")
		
	u.hyps <- as.integer(levels(src.uid)) # unique hypotheses
	nhyps <- length(u.hyps) # number of unique hypotheses

	results <- data.frame(matrix(0, 2*nhyps, 12), stringsAsFactors = F)
	colnames(results) = c('uid', 'name', 'regulation', 'correct.pred', 'incorrect.pred', 'score',
	'total.reachable', 'significant.reachable', 'total.ambiguous', 'significant.ambiguous',
	'unknown', 'ternary.pval')
	
	if (verbose)
		cat("\n Total number of hypothesis to consider:", nhyps)
	
	results[,1] <- rep(u.hyps,each=2)
	results[,2] <- rep(ents$name[match(u.hyps,ents$uid)],each=2)
	results[,3] <- rep(c("up","down"), nhyps)
	
	tab <- mapply(table,child.sgn,child.val)
	dim(tab) <- c(3,3,nhyps)
	dimnames(tab) <- list(sgn = c(-1,0,1), val = c(-1,0,1), uhyp = u.hyps)
	npp <- tab["1","1",]
	npm <- tab["1","-1",]
	npz <- tab["1","0",]
	nmp <- tab["-1","1",]
	nmm <- tab["-1","-1",]
	nmz <- tab["-1","0",]
	nrp <- tab["0","1",]
	nrm <- tab["0","-1",]
	nrz <- tab["0","0",]
	
	tab <- sapply(non.child.val,table)
	rownames(tab) <- c(-1,0,1)
	nzp <- tab["1",]
	nzm <- tab["-1",]
	nzz <- tab["0",]
	pvalt <- mapply(runCRE, npp=npp, npm=npm, npz=npz, nmp=nmp, nmm=nmm, nmz=nmz, 
		nrp=nrp, nrm=nrm, nrz=nrz, nzp=nzp, nzm=nzm, nzz=nzz, method = 'Ternary', 
		SIMPLIFY=FALSE)	
	qPlus  <- npp + npm + npz
	qMinus <- nmp + nmm + nmz
	qR 	   <- nrp + nrm + nrz
	qZero  <- nzp + nzm + nzz
	
	odd  <- seq.int(1L, len=nhyps, by=2L)
	even <- seq.int(2L, len=nhyps, by=2L)
	
	results[odd,  4] <- npp + nmm
	results[even, 4] <- nmp + npm
	results[odd,  5] <-npm + nmp
	results[even, 5] <- npp + nmm
	results[odd,  6] <-npp + nmm - (npm + nmp)
	results[even, 6] <- nmp + npm - (npp + nmm)
	results[	   ,  7] <- rep(qPlus + qMinus + qR, each = 2)
	results[   ,	  8] <- rep(npp + npm + nmp + nmm + nrp + nrm, each = 2) 
	results[   ,  9] <- rep(qR, each = 2)
	results[   , 10] <- rep(nrp + nrm, each = 2)
	results[	   ,	 11] <- rep(qZero, each = 2)
	results[odd, 12] <- sapply(pvalt,"[[","pval.up")
	results[even,12] <- sapply(pvalt,"[[","pval.down")
	
	return(results)
	
}
