\name{cvSGL}
\alias{cvSGL}

\title{Fit and cross-validate a GLM with a combination of lasso and group lasso with overlap regularization}

\description{Fits and cross-validates a regularized generalized linear model via penalized maximum likelihood.  The model is fit for a path of values of the penalty parameter, and a parameter value is chosen by cross-validation. Fits linear and logistic models.}

\usage{
cvSGL(data, index = NULL, weights=NULL, type = c("linear","logit"), alphas = seq(0,1,.1),
	nlam = 20, standardize = c("train","self","all","no"), nfold = 10, measure = c("ll","auc"), 
	maxit = 1000, thresh = 0.001, min.frac = 0.05, gamma = 0.8, step = 1, reset = 10, ncores = 1,
	lambdas = NULL, verbose = FALSE)
}

\arguments{
  \item{data}{A list with components $x$, an input matrix of dimension $(n,p)$, and $y$, a response vector of length $n$. For \code{type="logit"} $y$ should be binary}
  \item{index}{A $p$-vector indicating group membership of each covariate}
  \item{weights}{Optional vector of weights for the group penalties}
  \item{type}{Model type: "linear" or "logit"}
  \item{alphas}{Vector of mixing parameters. \code{alpha} = 1 is the lasso penalty.}
  \item{nlam}{Number of lambda values in the regularization path}
  \item{standardize}{Type of standardization for full data and CV folds.}
  \item{nfold}{Number of folds of the cross-validation loop}
  \item{measure}{Performance measure used to select the best values \code{alphas} and \code{lambdas}}
  \item{maxit}{Maximum number of iterations to convergence}
  \item{thresh}{Convergence threshold for change in beta}
  \item{min.frac}{Minimum value of the penalty parameter, as a fraction of the maximum value}
 \item{gamma}{Fitting parameter used for tuning backtracking (between 0 and 1)}
\item{step}{Fitting parameter used for initial backtracking step size (between 0 and 1)}
 \item{reset}{Fitting parameter used for taking advantage of local strong convexity in Nesterov momentum (number of iterations before momentum term is reset)}
 \item{ncores}{Number of computer cores to use in computations}
\item{lambdas}{User-specified sequence of lambda values for fitting. We recommend leaving this NULL and letting cvSGL self-select values}
  \item{verbose}{Logical flag for whether or not step number will be output}
}

\details{
The function executes \code{SGL} \code{nfold}+1 times; the initial run is to find the \code{lambda} sequence, subsequent runs are used to compute the cross-validated error rate and its standard deviation. By default, \code{weights} are the square roots of group sizes.
}

\value{An object of class \code{"cv.creNet"} and \code{"creNet"} with components 
	
	\item{fit}{The fitted model using the best values of \code{alphas} and \code{lambdas} (class \code{"creNet"})}
	\item{best.lambda}{Index and value of the best element in \code{lambdas}}
	\item{best.alpha}{Index and value of the best element in \code{alphas}}
  \item{lldiff}{
Cross-validation (negative) log likelihood for all \code{alphas} and \code{lambdas} (=squared error loss if \code{type=linear})}
  \item{llSD}{Approximate standard deviations of \code{lldiff}}
  \item{AUC}{Area Under the Curve}
  \item{lambdas}{Values of \code{lambda} used in cross-validation.}
  \item{alphas}{User-specified argument \code{alphas}.}

}



\references{Simon, N., Friedman, J., Hastie, T., and Tibshirani, R. (2011)
  \emph{A Sparse-Group Lasso}, \cr
  \url{http://web.stanford.edu/~hastie/Papers/SGLpaper.pdf}
  }
\author{Kourosh Zarringhalam and David Degras\cr

Modified from SGL package: Noah Simon, Jerome Friedman, Trevor Hastie, and Rob Tibshirani

Maintainer: Kourosh Zarringhalam <kourosh.zarringhalam@umb.edu>
}
\seealso{\code{creSGL}
}
\examples{
set.seed(1)
n = 50; p = 100; size.groups = 10
index <- ceiling(1:p / size.groups)
X = matrix(rnorm(n * p), ncol = p, nrow = n)
beta = (-2:2)
y = X[,1:5] \%*\% beta + 0.1*rnorm(n)
data = list(x = X, y = y)
weights = rep(1, size.groups)
cvFit = cvcreSGL(data, index, weights, type = "linear", maxit = 1000, thresh = 0.001, min.frac = 0.05, nlam = 100, gamma = 0.8, nfold = 10, standardize = TRUE, verbose = FALSE, step = 1, reset = 10, alpha = 0.05, lambdas = NULL)

}
\keyword{model}
\keyword{regression}
