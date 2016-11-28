SGL <- function(data, index=NULL, weights=NULL, type = c("linear","logit"), alphas = 0.95, nlam = 20,
	standardize = TRUE, maxit = 1000, thresh = 0.001, min.frac = 0.1, gamma = 0.8, step = 1, reset = 10,
	ncores = 1, lambdas = NULL)
{
	
	X.transform <- NULL
	nalpha <- length(alphas)
	if (nalpha > 1 & !is.list(lambdas))
			lambdas <- rep_len(list(lambdas),nalpha)
	
	if (standardize == TRUE) {
		stand <- standardize(data$x)
		data$x <- stand$x
		X.transform <- list(X.scale = stand$scale, X.means = stand$center)		
	}
	
	if (is.null(index)) 
		index <- rep(1,ncol(data$x))
	stopifnot(length(index) == ncol(data$x))
	
	type <- match.arg(type)
	oneD <- switch(type, linear = oneDim, logit = oneDimLogit)
		
	if (is.null(weights)) {
		weights <- as.numeric(sqrt(table(index)))
	} else if (length(weights)==1) {
	weights <- rep_len(weights,length(unique(index)))
	}
		
	if (type == "linear") {
		intercept <- ifelse(standardize, mean(data$y), 0)
		data$y <- data$y - intercept 
	}

	if (nalpha == 1) {	
		Sol <- oneD(data, index, weights, thresh, inner.iter = maxit, 
			outer.iter = maxit, outer.thresh = thresh, min.frac = min.frac, 
			nlam = nlam, lambdas = lambdas, gamma = gamma, 
			step = step, reset = reset, alpha = alphas)

		Sol$alphas <- alphas 
		Sol$type <- type
		if (type == "linear") Sol$intercept <- intercept
		Sol$X.transform <- X.transform
		class(Sol) <- "creNet"

	} else { # multiple alpha values
		Sol <- mcmapply(oneD, alpha = as.list(alphas), lambdas = lambdas, 
				MoreArgs = list(data=data, index=index, weights=weights, thres=thresh, inner.iter=maxit, 
				outer.iter=maxit, outer.thresh=thresh, min.frac = min.frac, nlam=nlam, gamma=gamma, 
				step=step, reset=reset), SIMPLIFY = FALSE, mc.cores = ncores)

		for (a in 1:nalpha) {
			Sol[[a]]$type <- type
			Sol[[a]]$alpha <- alphas[a]
			if (type == "linear") Sol[[a]]$intercept <- intercept
			Sol[[a]]$X.transform <- X.transform
			class(Sol[[a]]) <- "creNet"
		}					
	}

	return(Sol)
}
	




		
cvSGL <- function(data, index = NULL, weights=NULL, type = c("linear","logit"), alphas = seq(0,1,.1),
	nlam = 20, standardize = c("train","self","all","no"), nfold = 10, measure = c("ll","auc"), num.iter = 10, 
	maxit = 1000, thresh = 0.001, min.frac = 0.05, gamma = 0.8, step = 1, reset = 10, ncores = 1,
	lambdas = NULL, verbose = FALSE)
{
	
	type <- match.arg(type)
	standardize <- match.arg(standardize)
	X.transform <- NULL
	nalpha <- length(alphas)
	stand.in <- (standardize %in% c("train","self")) # flag for standardization of training folds
	stand.out <- ifelse(standardize == "self","self","train") # standardization method for holdout fold
	verbose.cv <- (verbose & ncores > 1) # limited printout if parallel computing
	
	
	best.alpha = best.lambda = matrix(0, nrow = num.iter, ncol = 2)
	colnames(best.alpha) = colnames(best.lambda) = c("index","value")
	
	## AUC related functions
	f.auc  <- function(probs) roc(probs,data$y)$AUC
	f.dist <- function(probs) roc(probs,data$y)$opt.thresh.dist
	f.f1   <- function(probs) roc(probs,data$y)$opt.thresh.f1
	f.ba   <- function(probs) roc(probs,data$y)$opt.thresh.ba
	opt.thresh.dist = opt.thresh.f1 = opt.thresh.ba = numeric(num.iter)
	
	# if (measure == 'auc') {
	#   f.auc  <- function(probs) roc(probs,data$y)$AUC
	#   f.dist <- function(probs) roc(probs,data$y)$opt.thresh.dist
	#   f.f1   <- function(probs) roc(probs,data$y)$opt.thresh.f1
	#   f.ba   <- function(probs) roc(probs,data$y)$opt.thresh.ba
	#   opt.thresh.dist = opt.thresh.f1 = opt.thresh.ba = numeric(num.iter)
	# }
	# 
	## Standardized data 
	
	if (standardize == "no") {
		X <- data$x
		X.transform <- NULL
		intercept <- 0
		y <- data$y
	} else {	
		stand <- standardize(data$x)
		X <- stand$x
		X.transform <- list(X.scale = stand$scale, X.means = stand$center)	
		intercept <- ifelse(type == "linear",mean(data$y),0)
		y <- data$y - intercept
		if (standardize == "all") {
			data$x <- X
			data$y <- y
			stand.out <- 'no'
		}
	}
		
	## Set the weights to square roots of group sizes if unspecified
	if (is.null(weights)) {
		weights <- as.vector(sqrt(table(index)))
	} else if (length(weights)==1) {
		weights <- rep_len(weights,length(unique(index)))
	}
	
	## Determine the regularization parameters lambda if not specified
 	if (is.null(lambdas)) {
		lambdas <- vector("list",nalpha)
	    for (a in 1:nalpha)
			lambdas[[a]] <- pathCalc(list(x=X,y=y), index, weights, alpha = alphas[a], 
				min.frac = min.frac, nlam = nlam, type = type)
		if (nalpha == 1) lambdas <- unlist(lambdas)		
	} else {
	## Sort and recycle lambdas if provided as single vector	
		if (nalpha > 1 & !is.list(lambdas)) { 
		lambdas <- sort(lambdas,decreasing=TRUE)
		lambdas <- rep(list(lambdas),nalpha)
		}
	}

	## Create the cross-validation folds
	n <- length(data$y)
	size = as.integer(floor(n/nfold))
	ind.split = c(seq.int(0L,by=size,len=nfold),n)
	
	## Result objects
	lambdas.list <- if (is.list(lambdas)) lambdas else list(lambdas)
	lldiffFold <- lapply(lambdas.list, function(x) matrix(0,nfold,length(x))) 
	pred <- lapply(lambdas.list, function(x) matrix(0,n,length(x))) # predicted responses 		
  pred <- lapply(1:num.iter, function(x) pred) ## repeate the list for all itterations

  AUCs <- lapply(lambdas.list, function(x) matrix(0,length(x),1))
  AUCs <- lapply(1:num.iter, function(x) AUCs)
  lldiffs <- lapply(lambdas.list, function(x) matrix(0,length(x),1))
  lldiffs <- lapply(1:num.iter, function(x) lldiffs)
  llSDs <- lapply(lambdas.list, function(x) matrix(0,length(x),1))
  llSDs <- lapply(1:num.iter, function(x) llSDs)

	##@@  CROSS-VALIDATION LOOP @@##
	for(iter in 1:num.iter){
	  if(verbose){
	    cat('\n\n')
	    cat(paste('running iteration', iter))
	    cat('\n')
	  }
	  
	  samp.counter = 1
	  while(TRUE & samp.counter < 100){
	    Flag = TRUE
	    ind <- sample(1:n, replace = F)
	    for (i in 1:nfold)
	    {
	      ind.out <- ind[(ind.split[i]+1):ind.split[i+1]] # holdout fold
	      ind.in <- setdiff(ind,ind.out) # training folds
	      if(sum(data$y[ind.in]==0) <= 3 || sum(data$y[ind.in]==1) <= 3){
	        Flag = FALSE
	        break
	      }
	    }
	    
	    if(Flag)
	      break
	    samp.counter = samp.counter + 1
	  }
	  
	  if(samp.counter >= 100){
	    cat('\n WARNING \n')
	    cat(' Numer of samples in of of the classes <= 3')
	  }
	    
	    
	  #ind <- sample(1:n, replace = FALSE)
	  
	  for (i in 1:nfold)
	  {
	    if (verbose)
	      cat("\n*** NFOLD ", i, "***")
	    ind.out <- ind[(ind.split[i]+1):ind.split[i+1]] # holdout fold
	    ind.in <- setdiff(ind,ind.out) # training folds
	    new.data <- list(x = data$x[ind.in,], y = data$y[ind.in]) # training data
	    
	    ## Fit model on training data 	
	    Sol <- SGL(data = new.data, index = index, weights = weights, type = type, maxit = maxit, 
	               thresh = thresh, min.frac = min.frac, lambdas = lambdas, standardize = stand.in,
	               gamma = gamma, step = step, reset = reset, alphas = alphas, ncores = ncores)
	    
	    
	    ## Predict responses for holdout fold
	    outer.x <- data$x[ind.out,,drop = F]
	    if (nalpha == 1) Sol <- list(Sol)
	    pred.y <- lapply(Sol, predict, newX=outer.x, lam=NULL, standardize=stand.out)[[1]]
	    for (a in 1:nalpha) {
	      pya <- pred.y[[a]]
	      
	      if (type == "linear") {
	        lldiffFold[[a]][i,] <- colSums((data$y[ind.out] - pya)^2)/2
	      } else { eps <- 1e-8
	      pya[pya < eps] <- eps
	      pya[pya > 1 - eps] <- 1 - eps
	      lldiffFold[[a]][i,] <- - colSums(data$y[ind.out] * log(pya) + 
	                                         (1-data$y[ind.out]) * log(1-pya)) 
	      } 
	      pred[[iter]][[a]][ind.out,] <- pya
	    }
	    
	  } # END cross-validation loop
	  
	  ## Compute cross-validation summaries
	  lldiff <- lapply(lldiffFold,colSums)
	  llSD <- lapply(lldiffFold, function(x) apply(x,2,sd) * sqrt(nfold))
	  
	  ## Find best parameters alpha and lambda	
	  AUC <- lapply(pred[[iter]], function(x) apply(x, 2, f.auc))
	  thresh.dist <- lapply(pred[[iter]], function(x) apply(x, 2, f.dist))
	  thresh.f1   <- lapply(pred[[iter]], function(x) apply(x, 2, f.f1))
	  thresh.ba   <- lapply(pred[[iter]], function(x) apply(x, 2, f.ba))
	  
	  if (measure == 'auc') {
	    ## There may be more than one max AUC. This can happen if
	    ## predicted probabilities are all close to 0 and 1 in which
	    ## case the AUC will be 0. In this case we take the largest lambda
	    best.alpha[iter, "index"] <- which.max(sapply(AUC,max)) ## For Alpha it is OK to pick smallest
	    tmp.AUC     <- AUC[[best.alpha[iter, "index"]]]
	    m.AUC       <- max(tmp.AUC)
	    m.AUC.ind   <- which(tmp.AUC == m.AUC)
	    opt.AUC.ind <- m.AUC.ind[length(m.AUC.ind)]
	    best.lambda[iter, "index"] <- opt.AUC.ind ## For lambda pick the largest
	    ##best.alpha[iter, "index"] <- which.max(sapply(AUC,max)) 
	    ##best.lambda[iter, "index"] <- which.max(AUC[[best.alpha[iter, "index"]]])
	  } else {
	    best.alpha[iter, "index"] <- which.min(sapply(lldiff,min)) 
	    best.lambda[iter, "index"] <- which.min(lldiff[[best.alpha[iter, "index"]]])
	  }
	  AUCs[[iter]] <- AUC
	  lldiffs[[iter]] <- lldiff
	  llSD[[iter]] <- llSD
	  
	  best.alpha[iter, "value"] <- alphas[best.alpha[iter,"index"]]
	  best.lambda[iter, "value"] <- lambdas.list[[best.alpha[iter, "index"]]][best.lambda[iter, "index"]]
	  opt.thresh.dist[iter] <- thresh.dist[[best.alpha[iter, "index"]]][best.lambda[iter, "index"]]
	  opt.thresh.f1[iter] <- thresh.f1[[best.alpha[iter, "index"]]][best.lambda[iter, "index"]]
	  opt.thresh.ba[iter] <- thresh.ba[[best.alpha[iter, "index"]]][best.lambda[iter, "index"]]
	  ## Unlist performance measures if single alpha value
	  # if (nalpha == 1) {
	  #   AUCs[[iter]] <- unlist(AUC)
	  #   lldiffs[[iter]] <- unlist(lldiff)
	  #   llSD[[iter]] <- unlist(llSD)
	  # }
	}
	
	## Re-run SGL with all data using best parameters lambda and alpha
	## (For best alpha, the full path of lambda is used to increase speed with warm starts)

	if (verbose == TRUE)
	  cat("\n*** Fit SGL on full dataset ***")
  
  Sols = list()
  for(iter in 1:num.iter){
    Sol <- SGL(data = list(x=X,y=y), index = index, weights = weights, type = type, maxit = maxit, 
               thresh = thresh, min.frac = min.frac, lambdas = lambdas.list[[best.alpha[iter, 1]]], ncores = 1L,
               standardize = FALSE, gamma = gamma, step = step, reset = reset, alphas = best.alpha[iter,2])
    Sol$beta <- Sol$beta[,best.lambda[iter,1],drop=FALSE]
    Sol$lambdas <- Sol$lambdas[best.lambda[iter,1]]
    Sol$intercept <- switch(type, linear = intercept, logit = Sol$intercept[best.lambda[iter,1]])
    Sol$X.transform <- X.transform
    Sol$opt.thresh.dist <- opt.thresh.dist[iter]
    Sol$opt.thresh.f1   <- opt.thresh.f1[iter]
    Sol$opt.thresh.ba   <- opt.thresh.ba[iter]
    Sols = c(Sols, list(Sol))
  }
	
	## Only keep result for best lambda	
	Sol <- list(fit = Sols, best.lambda = best.lambda, best.alpha = best.alpha,
		lldiffs = lldiffs, llSDs = llSDs, AUCs = AUCs, lambdas = lambdas, alphas = alphas)
	
	class(Sol) = c("cv.creNet", "creNet")
	
	return(Sol)
}
