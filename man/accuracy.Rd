\name{accuracy}
\alias{accuracy}

\title{Accuracy Measures Fort Predicted Responses}

\description{
This function computes accuracy measures for predicted numerical or binary responses.   
}

\usage{
assessModel(pred, obs, type=c("linear","logit"), verbose=TRUE, prednames=NULL)
}

\arguments{
  \item{pred}{Vector or matrix of predicted responses. Each column corresponds to a different lambda value}

  \item{obs}{Vector or observed binary responses}
  \item{type}{Type of prediction model}
  \item{verbose}{logical, defaults to \code{TRUE}}
  \item{prednames}{Optional vector of names for the different sets of predictions}
}

\details{
For linear regression, the function simply returns the root mean square error (RMSE) of the predicted responses. For logistic regression, the function computes the negative log-likelihood using the original vector \code{pred} whose values should be in $(0,1)$ (soft classification). For the purpose of hard classification, the function also converts 
this vector to binary values and computes the following classification measures: sensitivity a.k.a. true positive rate (TP/(TP+FN)), specificity a.k.a. true negative rate (TN/(FP+TN)), accuracy (TP+TN)/(TP+TN+FP+FN), and precision a.k.a. positive predictive value (TP/(TP+FP)). Two cutoff values are considered for the binary conversion: 0.5 and the proportion of 1's in \code{obs}.}

\value{
If \code{type="linear"}, the RMSE of the predictions.\cr
If \code{type="logit"}, a list with components 
\item{prior}{proportion of 1's in \code{obs}}
\item{ll}{negative log-likelihoods}
\item{classif}{classification measures for the cutoff value 0.5}
\item{classif.adj}{classification measures for the adjusted cutoff value \code{prior}}
}

\author{
David Degras
}

\seealso{
\code{\link{creSGL}}, \code{\link{predictSGL}}, \code{\link{runPredictions}}
}

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
	step = 1, reset = 10, alpha = 0.05, lambdas = NULL)

pred = predict(Fit, X, 5)
prednames = substr(fit$lambdas,1,4)
accuracy(pred, y, type="logit", prednames=prednames)
}

