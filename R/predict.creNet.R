predict.creNet <- function(obj, newX, lam=NULL, standardize = c("train","self", "no"))
{
  
	if ("cv.creNet" %in% class(obj)) obj <- obj$fit
	X <- newX
	if (!is.matrix(X)) X <- as.matrix(X)
	
	if(class(obj) != 'list') obj <- list(obj)
	pred.iter = list()
	for(iter in 1:length(obj)){
	  if (is.null(lam)) lam <- seq_along(obj[[iter]]$lambdas)
	  
	  ## Dimension checks
	  p <- NROW(obj[[iter]]$beta)
	  stopifnot(ncol(X) == p)
	  stopifnot(all(lam %in% seq_along(obj[[iter]]$lambdas)))
	  
	  ## Standardization
	  #standardize <- match.arg(standardize)
	  if (standardize == "train") {
	    X <- t(X) - obj[[iter]]$X.transform$X.means # transpose, center
	    X <- t(X / obj[[iter]]$X.transform$X.scale) # scale, transpose back
	  } else if (standardize == "self") {
	    X <- standardize(X)$x
	    # means <- colMeans(X)
	    # X <- t(X) - means # transpose, center
	    # vars <- sqrt(rowSums(X^2))
	    # X <- t(X / vars)  # scale, transpose back
	  }
	  
	  ## Predictor
	  intercept <- switch(obj[[iter]]$type, linear = obj[[iter]]$intercept, logit = obj[[iter]]$intercept[lam])
	  eta <- X %*% obj[[iter]]$beta[,lam] + intercept
	  
	  pred <- switch(obj[[iter]]$type, linear = eta, logit = exp(eta)/(1+exp(eta)))
	  pred.iter = c(pred.iter, list(pred))
	}

	return(pred.iter)
}
