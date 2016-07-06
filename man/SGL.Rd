\name{SGL}
\alias{SGL}

\title{Fit a GLM with a combination of lasso and group lasso with overlap regularization}
\description{Fit a regularized generalized linear model via penalized maximum likelihood.  The model is fit for a path of values of the penalty parameter. Fits  linear and logistic models.}
\usage{
SGL(data, index=NULL, weights=NULL, type = c("linear","logit"), alphas = 0.95, nlam = 20,
	standardize = TRUE, maxit = 1000, thresh = 0.001, min.frac = 0.1, gamma = 0.8, step = 1, reset = 10,
	ncores = 1, lambdas = NULL)
}

\arguments{
  \item{data}{A list with components $x$, an input matrix of dimension $(n,p)$, and $y$, a response vector of length $n$. For \code{type="logit"} $y$ should be binary.
}
  \item{index}{A $p$-vector indicating group membership of each covariate}
  \item{weights}{Optional vector of weights for the group penalties}
  \item{type}{Model type: "linear" or "logit"}
\item{alphas}{Vector of mixing parameters. The value 1 yields the lasso penalty. The value 0 yields the group lasso penalty.}
  \item{nlam}{Number of lambda values to use in the regularization path}
  \item{standardize}{Logical flag for variable standardization prior to fitting the model.}
  \item{maxit}{Maximum number of iterations to convergence}
  \item{thresh}{Convergence threshold for change in beta}
  \item{min.frac}{Minimum value of the penalty parameter, as a fraction of the maximum value}
 \item{gamma}{Fitting parameter used for tuning backtracking (between 0 and 1)}
\item{step}{Fitting parameter used for initial backtracking step size (between 0 and 1)}
 \item{reset}{Fitting parameter used for taking advantage of local strong convexity in Nesterov momentum (number of iterations before momentum term is reset)}
\item{lambdas}{Optional sequence of lambda values for fitting. We recommend leaving this NULL and letting SGL self-select values}
}

\details{
The sequence of models along the regularization path is fit by accelerated generalized gradient descent. If specified, argument \code{weights} should be a vector whose length is the number of groups. By default, \code{weights} is set to the square root of group sizes. 
}


\value{A single object of class \code{"creNet"} or a list of such objects if \code{alphas} has length >1.
	For each value in \code{alphas}, the result has components:
\item{beta}{A p by \code{nlam} matrix, giving the penalized MLEs for the nlam different models, where the index corresponds to the penalty parameter \code{lambda}}
\item{lambdas}{The actual sequence of \code{lambda} values used (penalty parameter)}
\item{type}{Response type (logistic/linear)}
\item{intercept}{For some model types, an intercept is fit}
\item{X.transform}{A list used in \code{predict} which gives the empirical mean and variance of the x matrix used to build the model}
\item{lambdas}{The sequence of lambda values used for fitting}
}

\references{Simon, N., Friedman, J., Hastie, T., and Tibshirani, R. (2011)
  \emph{A Sparse-Group Lasso}, \cr
  \url{http://web.stanford.edu/~hastie/Papers/SGLpaper.pdf}
}

\author{Kourosh Zarringhalam and David Degras\cr

Modified from SGL package: Noah Simon, Jerome Friedman, Trevor Hastie, and Rob Tibshirani

Maintainer: Kourosh Zarringhalam <kourosh.zarringhalam@umb.edu>
}

\seealso{\code{\link{cvSGL}}, \code{\link{predict.SGL}}}

\examples{
n = 50; p = 100; size.groups = 10
index <- ceiling(1:p / size.groups)
X = matrix(rnorm(n * p), ncol = p, nrow = n)
beta = (-2:2)
y = X[,1:5] \%*\% beta + 0.1*rnorm(n)
y = ifelse((exp(y) / (1 + exp(y))) > 0.5, 1, 0)
data = list(x = X, y = y)
weights = rep(1, size.groups)
fit = creSGL(data, index, weights, type = "linear", maxit = 1000, thresh = 0.001, 
	min.frac = 0.05, nlam = 100, gamma = 0.8, standardize = TRUE, verbose = FALSE, 
	step = 1, reset = 10, alphas = 0.05, lambdas = NULL)
}

\keyword{model}
\keyword{regression}
