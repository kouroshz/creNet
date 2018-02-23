pathCalc <- function(data, index, weights, alpha = 0.95, min.frac = 0.05,
                     nlam = 100, type = 'logit')
{
  X <- data$x
  y <- data$y
  n <- nrow(X)
  p <- ncol(X)
  ord <- order(index)
  index <- index[ord]
  X <- X[, ord]
  groups <- unique(index)
  num.groups <- length(groups)
  group.length <- as.vector(table(index))
  group.range <- c(0,cumsum(group.length))

  if (type == "linear") {
    resp <- y
  }
  if (type == "logit") {
    m.y <- mean(y)
    resp <- m.y * m.y * (1 - m.y) - (y - m.y)
  }

  solver <- function(x,alpha,w)
  {
    if (x[1] == 0 || w == 0) return(0)
    p <- length(x) + 1
    x[p] <- 0
    l <- (1:p) * x^2 - 2 * x * cumsum(x) + cumsum(x^2)
    r <- x^2 * ((1-alpha)*w/alpha)^2
    if (all(l <= r)) {
      return(0) } else i <- which.max(l > r) - 1L
    a <- i * alpha^2 - ((1-alpha)*w)^2
    if (a == 0) return(x[i]/alpha)
    b <- alpha * sum(x[1:i])
    c <- sum(x[1:i]^2)

    if(abs(b^2 - a * c) <= 1e-14)
      return(b/a)

    sol <- (b + c(-1,1) * sqrt(b^2 - a * c))/a
    valid <- which(alpha * sol <= x[i] & alpha * sol >= x[i+1])
    sol <- switch(length(valid)+1L, x[i]/alpha, sol[valid], sol[1])
    return(sol)
  }

  cors <- lapply(1:num.groups,
                 function(i) abs(crossprod(X[,(group.range[i]+1):group.range[i+1]], resp/n)))

  if ((alpha != 0) & (alpha != 1)) {
    cors <- lapply(cors, sort, decreasing=TRUE)
    lambda.max <- unlist(lapply(cors,'[[',1))/alpha
    active <- rep(TRUE,num.groups)
    max.lam <- 0
    for (i in 1:num.groups) {
      if (!active[i]) next
      lambda.max[i] <- solver(cors[[i]], alpha, weights[i])
      if (lambda.max[i] > max.lam) {
        max.lam <- lambda.max[i]
        active[lambda.max <= max.lam] <- FALSE
      }
    }
  }

  if (alpha == 1) {
    max.lam <- max(unlist(cors))
  }

  if (alpha == 0) {
    active <- which(weights > 0)
    if (!length(active)) return(rep(0,nlam))
    lambda.max <- rep(0,num.groups)
    lambda.max[active] <- sqrt(unlist(lapply(cors[active],crossprod))) /
      weights[active]
    max.lam <- max(lambda.max)
  }

  min.lam <- min.frac * max.lam
  lambdas <- exp(seq(log(max.lam), log(min.lam), len = nlam))
  return(lambdas)
}
