\name{predict.creNet}
\alias{predict.creNet}
\title{Predict Responses for New Observations Using Fitted creSGL Model}
\description{
Outputs predicted response values for new user input observations at a specified \code{lambda} values}

\usage{
predict(obj, newX, lam=NULL, standardize=c("train","self"))}

\arguments{
  \item{obj}{A fitted object of class \code{creSGL} or \code{cv.creSGL}}
  \item{newX}{A matrix of new covariates}
  \item{lam}{The indexes of lambda values to be used for prediction}
  \item{standardize}{The standardization method for \code{newX}}
  
}

\details{The matrix \code{newX} should have the same number of columns as the covariate matrix used to obtain \code{obj}. The argument \code{standardize} specifies how \code{newX} should be standardized. For the choice \code{"train"}, the means and variances attached to \code{obj} are used. For the choice \code{"self"}, the own means and variances of \code{newX} are used.}

\author{Kourosh Zarringhalam and David Degras\cr

Modified from SGL package: Noah Simon, Jerome Friedman, Trevor Hastie, and Rob Tibshirani

Maintainer: Kourosh Zarringhalam <kourosh.zarringhalam@umb.edu>
}


\references{Simon, N., Friedman, J., Hastie T., and Tibshirani, R. (2011)
  \emph{A Sparse-Group Lasso}, \cr
  \url{http://web.stanford.edu/~hastie/Papers/SGLpaper.pdf}
}

\value{A vector or matrix of predicted responses. Each column corresponds to a value of lambda.}

\seealso{\code{\link{SGL}} and \code{\link{cvSGL}}.}


\examples{
n = 50; p = 100; size.groups = 10
index <- ceiling(1:p / size.groups)
X = matrix(rnorm(n * p), ncol = p, nrow = n)
beta = (-2:2)
y = X[,1:5] \%*\% beta + 0.1*rnorm(n)
y = ifelse((exp(y) / (1 + exp(y))) > 0.5, 1, 0)
data = list(x = X, y = y)
weights = rep(1, size.groups)
Fit = creSGL(data, index, weights, type = "logit", maxit = 1000, thresh = 0.001, 
		min.frac = 0.05, nlam = 100, gamma = 0.8, standardize = TRUE, verbose = FALSE,
		step = 1, reset = 10, alphas = 0.05, lambdas = NULL)
X.new = matrix(rnorm(n * p), ncol = p, nrow = n)
predict(Fit, X.new, 5)
}

\keyword{models}
\keyword{regression}
