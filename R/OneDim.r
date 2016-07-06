oneDim <- function(data, index, weights, thresh = 0.0001, nlam = 20, 
	lambdas = NULL, beta.naught = rep(0,ncol(data$x)), inner.iter = 100, 
	outer.iter = 100, outer.thresh = 0.0001, gamma = 0.8, step = 1, reset = 10, 
	alpha = 0.95, min.frac = 0.05)
{

	if (is.null(lambdas)) {
		lambdas <- pathCalc(data = data, index = index, weights=weights,
			alpha=alpha, min.frac = min.frac, nlam = nlam, type = "linear")
	} else {
		nlam <- length(lambdas)
		lambdas <- sort(lambdas, decreasing = TRUE)
	}
	
	# X <- data$x
	# y <- data$y
	n <- nrow(data$x)
	p <- ncol(data$x)
	
	## Setting up group lasso stuff ##
	
	ord <- order(index)
	index <- index[ord]
	# X <- X[,ord]
	data$x <- data$x[,ord]
	dim(data$x) <- NULL
	unOrd <- match(1:length(ord),ord)
	
	## Coming up with other C++ info ##
	
  groups <- unique(index)
  num.groups <- length(groups)
  group.length <- as.vector(table(index))
  range.group.ind <- c(0,cumsum(group.length))

	
	## DONE SETTING UP C STUFF ##
	
	#alpha <- sqrt(2*log(p))/(1+sqrt(2*log(num.groups)/min(group.length)) + sqrt(2*log(p)))
	
	beta.old <- rep(0,p)
	beta.is.zero <- rep(1,num.groups)
	beta <- array(0, c(p,nlam))
	eta <- rep(0,n)
	
	for(k in 1:nlam) {
	
		## Commented out for warm start
		# beta.is.zero <- rep(1,num.groups)
		# beta.old <- rep(0,p)
		# eta <- rep(0,n)
	
		junk <- .C("linNest", X = as.double(data$x), y = as.double(data$y), 
			index = as.integer(index), nrow = as.integer(n), 
			ncol = as.integer(p), numGroup = as.integer(num.groups), 
			rangeGroupInd = as.integer(range.group.ind), 
			groupLen = as.integer(group.length), weights = as.double(weights), 
			lambda1 = as.double(lambdas[k]*alpha), 
			lambda2 = as.double(lambdas[k]*(1-alpha)), beta = as.double(beta.old),
			innerIter = as.integer(inner.iter), outerIter = as.integer(outer.iter), 
			thresh = as.double(thresh), outerThresh = as.double(outer.thresh), 
			eta = as.double(eta), gamma = as.double(gamma), 
			betaIsZero = as.integer(beta.is.zero), step = as.double(step), 
			reset = as.integer(reset))
	 
		beta.new <- junk$beta
		beta[,k] <- beta.new
		beta.is.zero <- junk$betaIsZero
		eta <- junk$eta
		beta.old <- beta.new
		
	}
	
	return(list(beta = beta[unOrd,], lambdas = lambdas))
	}
