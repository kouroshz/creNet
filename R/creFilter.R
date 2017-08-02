creFilter <- function(ents, rels, x.train, y.train, x.test = NULL, y.test = NULL, RNAseq = FALSE,
                      cre.sig = 0.01, de.sig = 0.01, type.weight = c("cre","log.cre","sqrt","both","none"), 
                      filter = TRUE, verbose = TRUE)
{
  
  type.weight <- match.arg(type.weight)
  
  ## Dimension checks
  stopifnot(NROW(x.train) == length(y.train))
  stopifnot(is.null(x.test) || NROW(x.test) == length(y.test))
  
  gene.ids <- as.numeric(colnames(x.train))
  pval <- weights <- NULL
  if (verbose)
    cat("\n Dimension of the data:\n", dim(x.train))
  
  L <- getRegulatorSliceIndex(ents, rels, gene.ids)
  u.hyps <- L$u.hyps
  child.uid <- L$child.uid
  child.sgn <- L$child.sgn
  groups <- L$groups
  child.id <- L$child.id
  data.slice <- L$data.slice
  
  
  ## Case where CRE is not required
  if (!filter && type.weight %in% c("sqrt","none")) {
    nhyps <- length(u.hyps)
    sig.ind <- 1:nhyps
    cre.priors.sig <- NULL
    sig.hyps <- NULL
    no.weights <- rep(1, nhyps)
    sqrt.weights <- sqrt(rle(groups)$lengths)
    ##sqrt.weights <- sqrt(sapply(groups,length))
    cre.weights <- NULL
    log.cre.weights <- NULL		
    both.weights <- NULL
  } 
  
  ## Other cases: CRE
  if (filter || type.weight %in% c("log.cre","cre","both")) {
    
    if (verbose) {
      cat("\n===================================\n")
      cat("\nRunning CRE")
    }
    ## Get differentially expressed genes
    L <- getDEGs(ents, x.train, y.train, de.sig, RNAseq)
    evidence <- L$evidence
    ## Get CRE priors
    cre.priors <- generatePriors(ents, rels, evidence, verbose)
    
    ## For each regulator, retain row ('up' or 'down')
    ## with smallest p-value and discard the other
    nhyps <- nrow(cre.priors)/2
    #pval <- fdrtool(cre.priors$ternary.pval, statistic = 'pvalue', plot = F, verbose = F)$qval
    pval <- cre.priors$ternary.pval
    odd <- seq.int(1L,by=2L,len=nhyps)
    even <- seq.int(2L,by=2L,len=nhyps)
    ind <- ifelse(pval[even] < pval[odd], even, odd)
    cre.priors <- cre.priors[ind,]
    pval <- pval[ind]
    
    ## Look at significant hypotheses only
    ind <- which(pval < cre.sig)
    if (length(ind) == 0 | !filter) ind <- 1:nhyps
    cre.priors.sig <- cre.priors[ind,]
    ## Get the significant ones
    sig.hyps <- cre.priors.sig$uid
    sig.ind <- match(sig.hyps,u.hyps)
    
    if (verbose) {
      cat("\n Number of significant hypotheses:", length(sig.hyps))
      cat("\n===================================")
      cat("\n Significant hypotheses", sig.hyps)
      cat("\n===================================")
      cat("\n Only considering genes down stream of the significant hypotheses")
    }
    
    ## Weights		
    cre.weights <- pval[sig.ind]/sum(pval[sig.ind])
    log.cre.weights <- 1/abs(log(pval[sig.ind]))		
    sqrt.weights <- sqrt(rle(groups)$lengths[sig.ind])
    ##sqrt.weights <- sqrt(sapply(groups[sig.ind],length))
    both.weights <- cre.weights * sqrt.weights
    no.weights <- rep(1,length(sig.hyps))
  }
  
  ## update child uid, id and sgn according to sig.hypes
  child.id = child.id[sig.ind]
  child.uid = child.uid[sig.ind]
  child.sgn = child.sgn[sig.ind]
  
  ## Duplicate data according to slices
  ## Note that slice is ordered according to u.hyps
  slice.ind <- unlist(data.slice[sig.ind])
  slice.train <- x.train[, slice.ind, drop=FALSE]
  x.train <- x.train[, unique(slice.ind), drop=FALSE]
  if (!is.null(x.test)) {
    slice.test <- x.test[, slice.ind, drop=FALSE]
    x.test <- x.test[, unique(slice.ind), drop=FALSE]
  } else slice.test <- NULL
  
  
  group.length <- sapply(data.slice[sig.ind],length)
  groups <- rep(1:length(sig.ind),group.length)	
  uid.groups <- rep(u.hyps[sig.ind], group.length)
  
  if (verbose) {
    cat("\n Dimension of the training data after filtering:\n", dim(slice.train))
    cat("\n Dimension of the training data (unique columns):\n", dim(x.train))
  }
  
  
  weights <- switch(type.weight, cre = cre.weights, log.cre = log.cre.weights, sqrt = sqrt.weights, 
                    both = both.weights, none = no.weights)
  
  L <- list(slice.train = slice.train, slice.test = slice.test, x.train = x.train, 
            x.test = x.test, u.hyps = u.hyps, group.length = group.length, sig.hyps = sig.hyps, 
            weights = weights, groups = groups, slice.ind = slice.ind, uid.groups = uid.groups, 
            cre.priors.sig = cre.priors.sig, child.uid = child.uid, child.sgn = child.sgn,
            cre.weights, log.cre.weights, sqrt.weights, both.weights)
  
  return(L)
  
}

creFilterOrig <- function(ents, rels, x.train, y.train, x.test = NULL, y.test = NULL, cre.sig = 0.01, de.sig = 0.01,
                          type.weight = 'cre', filter = TRUE, verbose = TRUE)
{
  
  stopifnot(NROW(x.train) == length(y.train))
  stopifnot(is.null(x.train) || NROW(x.test) == length(y.test))
  gene.ids <- as.numeric(colnames(x.train))
  pval <- weights <- NULL
  if (verbose == TRUE)
    cat("\n Dimension of the data:\n", dim(x.train))
  
  ## Get data slices, i.e, index of columns corresponding to each regulator
  ents.mRNA = ents[which(ents$type == 'mRNA'), ]
  
  ## For each hypothesis, identify the children and non-children indexes
  u.hyps = unique(rels$srcuid)
  child.uid = lapply(u.hyps, function(x) rels$trguid[which(rels$srcuid == x)])
  child.sgn = lapply(u.hyps, function(x) ifelse(rels$type[which(rels$srcuid == x)] == 'increase',
                                                1, ifelse(rels$type[which(rels$srcuid == x)] == 'decrease', -1, 0)))
  
  ## Get the data slices corresponding to each hypothesis
  child.id = lapply(child.uid, function(x) as.numeric(ents.mRNA$id[match(x,ents.mRNA$uid)])) ## to get the id
  data.slice = lapply(child.id, function(x) match(x, gene.ids))
  
  
  if(!filter & type.weight == 'none'){
    nhyps          <- length(u.hyps)
    sig.ind        <- 1:nhyps
    cre.priors.sig <- NULL
    sig.hyps       <- NULL
    none.weight    <- rep(1, nhyps)
    
  }else{
    ## CRE
    if (verbose == TRUE) {
      cat("\n===================================\n")
      cat("\nRunning CRE")
    }
    ## Get differentially expressed genes
    L <- getDEGs(ents, rels, x.train, y.train, de.sig, verbose)
    evidence <- L$evidence
    ## Get CRE priors
    cre.priors = generatePriors(ents, rels, evidence, verbose)
    ## For each regulator, retain row ('up' or 'down')
    ## with smallest p-value and discard the other
    nhyps <- nrow(cre.priors)/2
    pval <- fdrtool(cre.priors$ternary.pval,statistic = 'pvalue', plot = F, verbose = F)$qval
    ind <- numeric(nhyps)
    for (i in 1:nhyps)
      ind[i] <- ifelse(pval[2*i] < pval[2*i-1], 2*i, 2*i-1)
    cre.priors <- cre.priors[ind,]
    pval <- pval[ind]
    
    ## Look at significant hypothesis only
    ##type.weight <- 'cre'
    ind <- which(cre.priors$ternary.pval < cre.sig)
    
    if (length(ind) == 0 | !filter){
      ind <- 1:nhyps
    }
    
    cre.priors.sig = cre.priors[ind,]
    ## Get the significant ones
    sig.hyps = cre.priors.sig$uid
    if (verbose == TRUE) {
      cat("\n Number of significant hypotheses:", length(sig.hyps))
      cat("\n===================================")
      cat("\n Significant hypotheses", sig.hyps)
      cat("\n===================================")
      cat("\n Only considering genes down stream of the significant hypotheses")
    }
    ## Get the slices for significant hypothesis
    ## Note that slice is ordered according to u.hyps
    
    sig.ind <- which(u.hyps %in% sig.hyps)
    
    sig.u.hyp.ind = unlist(lapply(sig.hyps, function(x) which(u.hyps == x)))
    sig.u.hyp.child.id = lapply(sig.u.hyp.ind, function(x) child.id[x][[1]])
    sig.u.hyps.weights <- pval[sig.u.hyp.ind]/sum(pval[sig.u.hyp.ind])
    
    sqrt.weights <- unlist(lapply(child.id[sig.ind], function(x) sqrt(length(x))))
    both.weight <- sig.u.hyps.weights * sqrt.weights
    
  }
  
  
  ## Duplicate data according to slices
  ## Note that slice is ordered according to u.hyps
  slice.ind <- unlist(data.slice[sig.ind])
  slice.train <- x.train[,unlist(data.slice[sig.ind]),drop=FALSE]
  x.train <- x.train[,unique(unlist(data.slice[sig.ind])),drop=FALSE]
  if(!is.null(x.test)){
    slice.test <- x.test[,unlist(data.slice[sig.ind]),drop=FALSE]
    x.test <- x.test[,unique(unlist(data.slice[sig.ind])),drop=FALSE]
  }else{
    x.test <- NULL
    slice.test <- NULL
  }
  
  group.length <- sapply(data.slice[sig.ind], length)
  groups <- rep(1:length(u.hyps[sig.ind]), group.length)
  uid.groups <- rep(u.hyps[sig.ind], group.length)
  
  if (verbose == TRUE){
    cat("\n Dimension of the training data after filtering:\n", dim(slice.train))
    cat("\n Dimension of the training data (unique columns):\n", dim(x.train))
  }
  
  
  weights = switch(type.weight, cre = sig.u.hyps.weights, sqrt = sqrt.weights, both = both.weight, none = none.weight)
  
  L <- list(slice.train = slice.train, slice.test = slice.test, x.train = x.train, x.test = x.test, u.hyps = u.hyps,
            group.length = group.length,sig.hyps = sig.hyps, weights = weights, groups = groups, slice.ind = slice.ind,
            uid.groups = uid.groups, cre.priors.sig = cre.priors.sig)
  return(L)
  
}
# 
# 
# creFilter <- function(ents, rels, x.train, y.train, x.test = NULL, y.test = NULL, 
#                       cre.sig = 0.01, de.sig = 0.01, type.weight = c("cre","log.cre","sqrt","both","none"), 
#                       filter = TRUE, verbose = TRUE)
# {
#   
#   type.weight <- match.arg(type.weight)
#   
#   ## Dimension checks
#   stopifnot(NROW(x.train) == length(y.train))
#   stopifnot(is.null(x.test) || NROW(x.test) == length(y.test))
#   
#   gene.ids <- as.numeric(colnames(x.train))
#   pval <- weights <- NULL
#   if (verbose)
#     cat("\n Dimension of the data:\n", dim(x.train))
#   
#   L <- getRegulatorSliceIndex(ents, rels, gene.ids)
#   u.hyps <- L$u.hyps
#   child.uid <- L$child.uid
#   child.sgn <- L$child.sgn
#   groups <- L$groups
#   child.id <- L$child.id
#   data.slice <- L$data.slice
#   
#   
#   ## Case where CRE is not required
#   if (!filter && type.weight %in% c("sqrt","none")) {
#     nhyps <- length(u.hyps)
#     sig.ind <- 1:nhyps
#     cre.priors.sig <- NULL
#     sig.hyps <- NULL
#     no.weights <- rep(1, nhyps)
#     sqrt.weights <- sqrt(rle(groups)$lengths)
#     cre.weights <- NULL
#     log.cre.weights <- NULL		
#     sqrt.weights <- NULL
#     both.weights <- NULL
#   } 
#   
#   ## Other cases: CRE
#   if (filter | type.weight %in% c("log.cre","cre","both")) {
#     
#     if (verbose) {
#       cat("\n===================================\n")
#       cat("\nRunning CRE")
#     }
#     ## Get differentially expressed genes
#     L <- getDEGs(ents, x.train, y.train, de.sig)
#     evidence <- L$evidence
#     ## Get CRE priors
#     cre.priors <- generatePriors(ents, rels, evidence, verbose)
#     
#     ## For each regulator, retain row ('up' or 'down')
#     ## with smallest p-value and discard the other
#     nhyps <- nrow(cre.priors)/2
#     #pval <- fdrtool(cre.priors$ternary.pval, statistic = 'pvalue', plot = F, verbose = F)$qval
#     #score <- cre.priors$score
#     pval <- cre.priors$ternary.pval
#     odd <- seq.int(1L,by=2L,len=nhyps)
#     even <- seq.int(2L,by=2L,len=nhyps)
#     #ind <- ifelse(score[even] > score[odd], even, odd)
#     ind <- ifelse(pval[even] < pval[odd], even, odd)
#     cre.priors <- cre.priors[ind,]
# 
#     #score <- score[ind]
#     #pval <- cre.priors$ternary.pval
#     pval <- pval[ind]
#     
#     ## Look at significant hypotheses only
#     ind <- which(pval < cre.sig)
#     if (length(ind) == 0 | !filter) ind <- 1:nhyps
#     cre.priors.sig <- cre.priors[ind,]
#     ## Get the significant ones
#     sig.hyps <- cre.priors.sig$uid
#     sig.ind <- match(sig.hyps,u.hyps)
#     
#     if (verbose) {
#       cat("\n Number of significant hypotheses:", length(sig.hyps))
#       cat("\n===================================")
#       cat("\n Significant hypotheses", sig.hyps)
#       cat("\n===================================")
#       cat("\n Only considering genes down stream of the significant hypotheses")
#     }
#     
#     ## Weights		
#     cre.weights <- pval[sig.ind]/sum(pval[sig.ind])
#     log.cre.weights <- 1/abs(log(pval[sig.ind]))		
#     #sqrt.weights <- sqrt(sapply(groups[sig.ind],length))
#     sqrt.weights <- sqrt(rle(groups[sig.ind])$lengths)
#     both.weights <- cre.weights * sqrt.weights
#     no.weights <- rep(1,length(sig.hyps))
#     
#   }
#   
#   
#   ## Duplicate data according to slices
#   ## Note that slice is ordered according to u.hyps
#   slice.ind <- unlist(data.slice[sig.ind])
#   slice.train <- x.train[, slice.ind, drop=FALSE]
#   x.train <- x.train[, unique(slice.ind), drop=FALSE]
#   if (!is.null(x.test)) {
#     slice.test <- x.test[, slice.ind, drop=FALSE]
#     x.test <- x.test[, unique(slice.ind), drop=FALSE]
#   } else slice.test <- NULL
#   
#   
#   group.length <- sapply(data.slice[sig.ind],length)
#   groups <- rep(1:length(sig.ind),group.length)	
#   uid.groups <- rep(u.hyps[sig.ind], group.length)
#   child.uid <- child.uid[sig.ind]
#   child.sgn <- child.sgn[sig.ind]
#   if (verbose) {
#     cat("\n Dimension of the training data after filtering:\n", dim(slice.train))
#     cat("\n Dimension of the training data (unique columns):\n", dim(x.train))
#   }
#   
#   
#   weights <- switch(type.weight, cre = cre.weights, log.cre = log.cre.weights, sqrt = sqrt.weights, 
#                     both = both.weights, none = no.weights)
#   
#   L <- list(slice.train = slice.train, slice.test = slice.test, x.train = x.train, 
#             x.test = x.test, u.hyps = u.hyps, group.length = group.length, sig.hyps = sig.hyps, 
#             weights = weights, groups = groups, slice.ind = slice.ind, uid.groups = uid.groups, 
#             cre.priors.sig = cre.priors.sig,child.uid = child.uid,child.sgn = child.sgn,
#             cre.weights, log.cre.weights, sqrt.weights, both.weights)
#   
#   return(L)
#   
# }
# 
# creFilterOrig <- function(ents, rels, x.train, y.train, x.test = NULL, y.test = NULL, cre.sig = 0.01, de.sig = 0.01,
#                       type.weight = 'cre', filter = TRUE, verbose = TRUE)
# {
#   
#   stopifnot(NROW(x.train) == length(y.train))
#   stopifnot(is.null(x.train) || NROW(x.test) == length(y.test))
#   gene.ids <- as.numeric(colnames(x.train))
#   pval <- weights <- NULL
#   if (verbose == TRUE)
#     cat("\n Dimension of the data:\n", dim(x.train))
#   
#   ## Get data slices, i.e, index of columns corresponding to each regulator
#   ents.mRNA = ents[which(ents$type == 'mRNA'), ]
#   
#   ## For each hypothesis, identify the children and non-children indexes
#   u.hyps = unique(rels$srcuid)
#   child.uid = lapply(u.hyps, function(x) rels$trguid[which(rels$srcuid == x)])
#   child.sgn = lapply(u.hyps, function(x) ifelse(rels$type[which(rels$srcuid == x)] == 'increase',
#                                                 1, ifelse(rels$type[which(rels$srcuid == x)] == 'decrease', -1, 0)))
#   
#   ## Get the data slices corresponding to each hypothesis
#   child.id = lapply(child.uid, function(x) as.numeric(ents.mRNA$id[match(x,ents.mRNA$uid)])) ## to get the id
#   data.slice = lapply(child.id, function(x) match(x, gene.ids))
#   
#   
#   if(!filter & type.weight == 'none'){
#     nhyps          <- length(u.hyps)
#     sig.ind        <- 1:nhyps
#     cre.priors.sig <- NULL
#     sig.hyps       <- NULL
#     none.weight    <- rep(1, nhyps)
#     
#   }else{
#     ## CRE
#     if (verbose == TRUE) {
#       cat("\n===================================\n")
#       cat("\nRunning CRE")
#     }
#     ## Get differentially expressed genes
#     L <- getDEGs(ents, rels, x.train, y.train, de.sig, verbose)
#     evidence <- L$evidence
#     ## Get CRE priors
#     cre.priors = generatePriors(ents, rels, evidence, verbose)
#     ## For each regulator, retain row ('up' or 'down')
#     ## with smallest p-value and discard the other
#     nhyps <- nrow(cre.priors)/2
#     pval <- fdrtool(cre.priors$ternary.pval,statistic = 'pvalue', plot = F, verbose = F)$qval
#     ind <- numeric(nhyps)
#     for (i in 1:nhyps)
#       ind[i] <- ifelse(pval[2*i] < pval[2*i-1], 2*i, 2*i-1)
#     cre.priors <- cre.priors[ind,]
#     pval <- pval[ind]
#     
#     ## Look at significant hypothesis only
#     ##type.weight <- 'cre'
#     ind <- which(cre.priors$ternary.pval < cre.sig)
#     
#     if (length(ind) == 0 | !filter){
#       ind <- 1:nhyps
#     }
#     
#     cre.priors.sig = cre.priors[ind,]
#     ## Get the significant ones
#     sig.hyps = cre.priors.sig$uid
#     if (verbose == TRUE) {
#       cat("\n Number of significant hypotheses:", length(sig.hyps))
#       cat("\n===================================")
#       cat("\n Significant hypotheses", sig.hyps)
#       cat("\n===================================")
#       cat("\n Only considering genes down stream of the significant hypotheses")
#     }
#     ## Get the slices for significant hypothesis
#     ## Note that slice is ordered according to u.hyps
#     
#     sig.ind <- which(u.hyps %in% sig.hyps)
#     
#     sig.u.hyp.ind = unlist(lapply(sig.hyps, function(x) which(u.hyps == x)))
#     sig.u.hyp.child.id = lapply(sig.u.hyp.ind, function(x) child.id[x][[1]])
#     sig.u.hyps.weights <- pval[sig.u.hyp.ind]/sum(pval[sig.u.hyp.ind])
#     
#     sqrt.weights <- unlist(lapply(child.id[sig.ind], function(x) sqrt(length(x))))
#     both.weight <- sig.u.hyps.weights * sqrt.weights
#     
#   }
#   
#   
#   ## Duplicate data according to slices
#   ## Note that slice is ordered according to u.hyps
#   slice.ind <- unlist(data.slice[sig.ind])
#   slice.train <- x.train[,unlist(data.slice[sig.ind]),drop=FALSE]
#   x.train <- x.train[,unique(unlist(data.slice[sig.ind])),drop=FALSE]
#   if(!is.null(x.test)){
#     slice.test <- x.test[,unlist(data.slice[sig.ind]),drop=FALSE]
#     x.test <- x.test[,unique(unlist(data.slice[sig.ind])),drop=FALSE]
#   }else{
#     x.test <- NULL
#     slice.test <- NULL
#   }
#   
#   group.length <- sapply(data.slice[sig.ind], length)
#   groups <- rep(1:length(u.hyps[sig.ind]), group.length)
#   uid.groups <- rep(u.hyps[sig.ind], group.length)
#   
#   if (verbose == TRUE){
#     cat("\n Dimension of the training data after filtering:\n", dim(slice.train))
#     cat("\n Dimension of the training data (unique columns):\n", dim(x.train))
#   }
#   
#   
#   weights = switch(type.weight, cre = sig.u.hyps.weights, sqrt = sqrt.weights, both = both.weight, none = none.weight)
#   
#   L <- list(slice.train = slice.train, slice.test = slice.test, x.train = x.train, x.test = x.test, u.hyps = u.hyps,
#             group.length = group.length,sig.hyps = sig.hyps, weights = weights, groups = groups, slice.ind = slice.ind,
#             uid.groups = uid.groups, cre.priors.sig = cre.priors.sig)
#   return(L)
#   
# }