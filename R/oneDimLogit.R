oneDimLogit <- function(data, index, weights, thresh, lambdas, inner.iter, outer.iter,
                        outer.thresh, gamma, step, reset, alpha, min.frac, nlam)
{

  if (is.null(lambdas)) {
    lambdas <- pathCalc(data = data, index = index, weights=weights, alpha=alpha,
                        min.frac = min.frac, nlam = nlam, type = "logit")
  } else {
  	nlam <- length(lambdas)
  	lambdas <- sort(lambdas, decreasing = TRUE)
  }
  
  # X <- data$x
  # y <- data$y
  n <- NROW(data$x)
  p <- NCOL(data$x)

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

  beta <- matrix(0, nrow = p, ncol = nlam)
  eta <- rep(0,n)
  intercept <- rep(log(sum(data$y)) - log(n-sum(data$y)), nlam)
  eta = eta + intercept[1]
  beta.is.zero <- rep(1, num.groups)
  beta.old <- rep(0,p)

  for (k in 1:nlam) {

    junk <- .C("logitNest", X = as.double(data$x), y = as.integer(data$y),
               index = as.integer(index), nrow = as.integer(n),
               ncol = as.integer(p), numGroup = as.integer(num.groups),
               rangeGroupInd = as.integer(range.group.ind), groupLen = as.integer(group.length),
               weights = as.double(weights), lambda1 = as.double(alpha*lambdas[k]),
               lambda2 = as.double((1-alpha)*lambdas[k]), beta = as.double(beta.old),
               innerIter = as.integer(inner.iter), outerIter = as.integer(outer.iter),
               thresh = as.double(thresh), outerThresh = as.double(outer.thresh),
               eta = as.double(eta), gamma = as.double(gamma),
               betaIsZero = as.integer(beta.is.zero),
               betaZero = as.double(intercept[k]), step = as.double(step))


    intercept[k] = junk$betaZero

    if (k < nlam){
      intercept[k+1] = intercept[k]
    }
    beta.new <- junk$beta
    beta[,k] <- beta.new
    beta.is.zero <- junk$betaIsZero
    eta <- junk$eta
    beta.old <- beta.new

  }

  return(list(beta = beta[unOrd,], lambdas = lambdas, intercept = intercept))
}
