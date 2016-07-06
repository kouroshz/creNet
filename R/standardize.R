standardize <- function(x) {
  if (!is.matrix(x)) dim(x) <- c(length(x),1)
  means <- colMeans(x)
  x <- t(x) - means
  vars <- sqrt(rowSums(x^2))
  x <- t(x/vars)
  return(list(x=x, center=means, scale=vars))
}
