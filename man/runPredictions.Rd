\name{runPredictions}
\alias{runPredictions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit SGL Model, Predict Responses, and Assess Prediction Accuracy}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
runPredictions(ents, rels, x.train, y.train, x.test = NULL, y.test = NULL, type = c("linear", "logit"), alpha = 0.95, nlam = 20, min.frac = 0.05, cv = TRUE, nfold = 10, cre = c("filter", "weight", "both", "none"), cre.sig = 0.01, standardize = c("all","self","train","no"), cores = 1, verbose = TRUE)
}

\arguments{
  \item{ents}{Entry data frame typically created by \code{\link{processMicroarray}}}
 
  \item{rels}{Relation data frame typically created by \code{\link{processMicroarray}}}
  
  \item{x.train}{Vector of responses for training data}

  \item{y.train}{Matrix of covariates for training data}

  \item{x.test}{Optional matrix of covariates for testing data}

  \item{y.test}{Optional vector of responses for testing data}

  \item{type}{Type of regression model: \code{linear} or \code{logit}}

  \item{alpha}{Tradeoff between lasso penalty and group lasso penalty. 
  	\code{alpha=1} is pure lasso, \code{alpha=0} is pure group lasso. 
  	Several values can be specified }

  \item{nlam}{Number of lambda values for the regularization path}

  \item{min.frac}{Smallest lambda value as a fraction of the largest}

  \item{cv}{logical flag: should the data be cross-validated? }

  \item{nfold}{Number of folds for cross-validation}

  \item{cre}{CRE method for filtering and/or computing group weights}

  \item{cre.sig}{significance level for CRE filtering}

  \item{standardize}{type of standardization}
 
  \item{cores}{Number of cores to be used in computations. \code{cores>1} will result in parallel computations}
 
  \item{verbose}{logical flag for verbosity level}
}

\details{
The possible values of \code{standardization} are: \code{"all"}: training and testing data are concatenated and then standardized, \code{"self"}: each data set (training and testing) is standardized separately, \code{"train"}: both training and testing data are standardized using the means and scale of the training data, \code{"no"}: no standardization. 
}

\value{
A list with components

\item{fit}{Fitted object(s) of class \code{creSGL} if \code{cv=FALSE} or \code{cvcreSGL} if \code{cv=TRUE}; one fitted object per value of \code{alpha}}

\item{alpha}{Input argument \code{alpha}}

\item{bestlam}{Best value(s) of $lambda$ for cross-validation score, \code{NULL} otherwise; \code{bestlam} has same length as \code{alpha}}

\item{pred}{Vector/matrix of predictions for training data and for testing data if specified; each column corresponds to a value of \code{alpha}}

\item{accuracy}{Accuracy measures in prediction}

\item{slice}{Duplicated matrix of covariates for training data and for testing data if available}
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
Kourosh Zarringhalam and David Degras}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%%
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
