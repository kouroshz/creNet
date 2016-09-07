nested.cvSGL <- function(ents, rels, x, y, type = c("linear","logit"), alphas = seq(0,1,.1),
	nlam = 20, standardize = c("train","self","all","no"), nfold = 10, measure = c("ll","auc"), 
	type.weight = c("cre","log.cre","sqrt","both","none"), filter = TRUE, num.iter = 10,
	maxit = 1000, thresh = 0.001, min.frac = 0.05, gamma = 0.8, step = 1, reset = 10, cre.sig = 0.01,
	de.sig = 0.01, ncores = 1, lambdas = NULL, verbose = FALSE)
{

  type.weight <- match.arg(type.weight)
	standardize <- match.arg(standardize)
	measure <- match.arg(measure)
	n <- length(y)
	pred <- matrix(NA, nrow = n, ncol = num.iter)
	best.lambdas <- matrix(0, nrow = (nfold+1), ncol = 2*num.iter)
	best.alphas <- matrix(0, nrow = (nfold+1), ncol = 2*num.iter)
	
	size = as.integer(floor(n/(nfold+1)))
	ind.split.outer = c(seq.int(0L,by=size,len=(nfold+1)),n)
	outer.indecies <- lapply(1:(nfold+1), function(x) list()) # predicted responses 		
	outer.indecies <- lapply(1:num.iter, function(x) outer.indecies) ## repeate the list for all itterations
	

	## Optimal threshold values and correponding response labels
	opt.thresh.dist <- matrix(-1, nrow = (nfold + 1), ncol = num.iter)
	opt.thresh.f1   <- matrix(-1, nrow = (nfold + 1), ncol = num.iter)
	opt.thresh.ba   <- matrix(-1, nrow = (nfold + 1), ncol = num.iter)
	labs.ep <- matrix(NA, nrow = n, ncol = num.iter)
	labs.di <- matrix(NA, nrow = n, ncol = num.iter)
	labs.f1 <- matrix(NA, nrow = n, ncol = num.iter)
	labs.ba <- matrix(NA, nrow = n, ncol = num.iter)
	sig.groups = list()
	## generating CV folds.
	nonzero.genes = {}
	nonzero.coeffs = {}
	for(iter in 1:num.iter){
	  if(verbose){
	    cat('\n\n')
	    cat(paste('running outer iteration', iter))
	    cat('\n')
	  }

	  fold.counter = 1
	  while(TRUE & fold.counter < 100){
	    Flag = TRUE
	    ind <- sample(1:n)
	    for (o in 1:(nfold+1))
	    {
	      ind.outer <- ind[(ind.split.outer[o]+1):ind.split.outer[o+1]]
	      ind.inner <- setdiff(ind,ind.outer)
	      if(sum(y.train[ind.inner]==0) <= 3 || sum(y.train[ind.inner]==1) <= 3){
	        Flag = FALSE
	        break
	      }
	    }
	    
	    if(Flag)
	      break
	    
	    fold.counter = fold.counter + 1
	  }
	  
	  if(fold.counter > 100){
	    cat('\n Sample size in one classes is too small \n')
	    stop()
	  }
	  
	  for (o in 1:(nfold+1))
	  {
	    if (verbose)
	      cat("\n*** OUTER NFOLD", o, "***")
	    ind.outer <- ind[(ind.split.outer[o]+1):ind.split.outer[o+1]]
	    ind.inner <- setdiff(ind,ind.outer)
	    outer.indecies[[iter]][[o]] = list(ind.outer)
	    if (verbose)
	      cat("\n*** Running CRE ***")		
	    L <- creFilter(ents, rels, x[ind.inner,], y[ind.inner], cre.sig = cre.sig, de.sig = de.sig, 
	                   type.weight = type.weight, verbose = FALSE)
	    slice.train <- L$slice.train
	    slice.ind <- L$slice.ind
	    slice.test <- x[ind.outer, slice.ind]
	    groups <- L$groups
	    weights <- L$weights
	    uid.groups <- L$uid.groups
	    child.uid <- L$child.uid
	    
	    # if (standardize == 'all') {
	    # D <- standardize(rbind(slice.train,slice.test))
	    # ntrain <- length(ind.inner)
	    # slice.train <- D$x[1:ntrain, ]
	    # X.transform <- list(X.means = D$center, X.scale = D$scale)
	    # }
	    
	    inner.data <- list(x = slice.train, y = y[ind.inner])
	    # stand.in <- standardize
	    # if (standardize == "all") stand.in <- "no"
	    
	    ## Cross-validation
	    if (verbose)
	      cat("\n*** Running Cross-Validation ***")
		  inner.fit <- NULL
		  loop.counter = 1
	    while(is.null(inner.fit) & loop.counter < 10){
	      try(inner.fit <- cvSGL(data=inner.data, index=groups, weights=weights, type=type, alphas=alphas,num.iter = 1,
	                             nlam=nlam, standardize=standardize, nfold=nfold, measure=measure, maxit=maxit, thresh=thresh,
	                             min.frac=min.frac, gamma=gamma, step=step, reset=reset, ncores=ncores, lambdas=NULL,
	                             verbose=verbose))
	      loop.counter = loop.counter + 1
	    }
		  
		  if(loop.counter >= 10){
		    cat('\n cannot fit the model \n')
		    cat('\n is data size too small? \n')
		    stop()
		  }
	    # inner.fit <- cvSGL(data=inner.data, index=groups, weights=weights, type=type, alphas=alphas,num.iter = 1,
	    #                    nlam=nlam, standardize=standardize, nfold=nfold, measure=measure, maxit=maxit, thresh=thresh,
	    #                    min.frac=min.frac, gamma=gamma, step=step, reset=reset, ncores=ncores, lambdas=NULL,
	    #                    verbose=verbose)
	    
	    if (standardize == "all") {
	      slice <- rbind(slice.train,slice.test)
	      inner.fit$X.transform <- list(center = colMeans(slice), scale = apply(slice,2,sd))
	    }
	    
	    best.lambdas[o, (2*iter-1):(2*iter)] <- inner.fit$best.lambda
	    best.alphas[o, (2*iter-1):(2*iter)] <- inner.fit$best.alpha
	    stand.out <- ifelse(standardize == "self","self","train")
	    ##pred[ind.outer] <- predict(inner.fit,slice.test,NULL,stand.out)
	    ##
	    opt.thresh.dist[o, iter] <- inner.fit$fit[[1]]$opt.thresh.dist
	    opt.thresh.f1[o, iter]   <- inner.fit$fit[[1]]$opt.thresh.f1
	    opt.thresh.ba[o, iter]   <- inner.fit$fit[[1]]$opt.thresh.ba
	    ##
	    
	    ##
	    nonzero.genes.tmp = which(inner.fit$fit[[1]]$beta != 0)
	    nonzero.coeffs.tmp = inner.fit$fit[[1]]$beta[nonzero.genes.tmp]
	    sig.groups = c(sig.groups, list(ents[which(ents$uid %in% unique(uid.groups[unique(nonzero.genes.tmp)])),]))
	    nonzero.genes = c(nonzero.genes,  unlist(child.uid)[nonzero.genes.tmp])
	    nonzero.coeffs = c(nonzero.coeffs, nonzero.coeffs.tmp)
	    ##
	    
	    pred[ind.outer, iter] <- predict(inner.fit,slice.test,NULL,'self')[[1]]
	    
	    labs.ep[ind.outer, iter] <- ifelse(pred[ind.outer, iter] > 0.5, 1, 0)
	    labs.di[ind.outer, iter] <- ifelse(pred[ind.outer, iter] > inner.fit$fit[[1]]$opt.thresh.dist, 1, 0)
	    labs.f1[ind.outer, iter] <- ifelse(pred[ind.outer, iter] > inner.fit$fit[[1]]$opt.thresh.f1, 1, 0)
	    labs.ba[ind.outer, iter] <- ifelse(pred[ind.outer, iter] > inner.fit$fit[[1]]$opt.thresh.ba, 1, 0)
	  }
	  
	}

	L <- list(sig.hyps = unique(unlist(lapply(sig.groups, function(x) x$name))), pred = pred, 
	          opt.thresh.dist = opt.thresh.dist, opt.thresh.f1 = opt.thresh.f1, 
	          opt.thresh.ba = opt.thresh.ba,best.alphas = best.alphas, best.lambdas = best.lambdas,
	          labs.ep = labs.ep, labs.di = labs.di, labs.fq = labs.f1, labs.ba = labs.ba, 
	          outer.indecies = outer.indecies, nonzero.genes = nonzero.genes, nonzero.coeffs = nonzero.coeffs)
	
	class(L) = c("nested.cv.creNet")
	
	return(L)

}
