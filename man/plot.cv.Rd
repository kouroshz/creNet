\name{plot.cv.SGL}
\alias{plot.cv.SGL}
\title{plots the cross-validated error curve produced by cv.SGL}
\description{
Plots the cross-validated error curve, and confidence bounds for each \code{lambda} in our regularization path.}
\usage{
\method{plot}{cv.SGL}(x, ...)
}
\arguments{
  \item{x}{fitted \code{"cv.SGL"} object}
  \item{...}{additional arguments to be passed to plot}
}
\details{
A cross validated deviance plot is produced. More regularized models are to the right (less regularized to the left)}
\references{Simon, N., Friedman, J., Hastie T., and Tibshirani, R. (2011)
  \emph{A Sparse-Group Lasso}, \cr
  \url{http://www-stat.stanford.edu/~nsimon/creSGL.pdf}
}
\author{Kourosh Zarringhalam\cr

Modified from SGL package: Noah Simon, Jerome Friedman, Trevor Hastie, and Rob Tibshirani

Maintainer: Kourosh Zarringhalam <kourosh.zarringhalam@umb.edu>
}
\seealso{\code{creSGL} and \code{cv.SGL}.}
\examples{
n = 50; p = 100; size.groups = 10
index <- ceiling(1:p / size.groups)
X = matrix(rnorm(n * p), ncol = p, nrow = n)
beta = (-2:2)
y = X[,1:5] \%*\% beta + 0.1*rnorm(n)
y = ifelse((exp(y) / (1 + exp(y))) > 0.5, 1, 0)
data = list(x = X, y = y)
weights = rep(1, size.groups)
cvFit = cvSGL(data, index, weights, type = "logit", maxit = 1000, thresh = 0.001, min.frac = 0.05, nlam = 100, gamma = 0.8, nfold = 10, standardize = TRUE, verbose = FALSE, step = 1, reset = 10, alpha = 0.05, lambdas = NULL)
plot(cvFit)
}
\keyword{models}
\keyword{regression}
