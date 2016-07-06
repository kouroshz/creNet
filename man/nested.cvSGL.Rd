\name{nested.cvSGL}
\alias{nested.cvSGL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Nested Cross-validation for Sparse Group Lasso}
\description{
Accurate evaluation of the predictive performance of SGL}

\usage{
nested.cvSGL(ents, rels, x, y, type = c("linear", "logit"), alphas = seq(0, 1, 0.1), nlam = 20, 
	standardize = c("train", "self", "all", "no"), nfold = 10, measure = c("ll", "auc"), maxit = 1000, 
	thresh = 0.001, min.frac = 0.05, gamma = 0.8, step = 1, reset = 10, cre.sig = 0.01, de.sig = 0.01, 
	ncores = 1, lambdas = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ents}{data frame of KB entries}
  \item{rels}{data frame of KB relations}
  \item{x}{matrix of predictors (gene expression levels)}
  \item{y}{response vector (numeric for linear regression, binary for logistic regression)}
  \item{type}{type of model: linear or logistic}
  \item{alphas}{Vector of mixing parameters for lasso/group lasso penalty (0 = pure group lasso, 1 = pure lasso)}
  \item{nlam}{number of lambda values in the regularization path}
  \item{standardize}{type of data standardization to be performed}
  \item{nfold}{number of folds in the inner cross-validation}
  \item{measure}{performance measure used to select best lambdas and alphas: log-likelihood or area under curve}
  \item{maxit}{Maximum number of iterations to convergence}
  \item{thresh}{Convergence threshold for change in beta}
  \item{min.frac}{Minimum value of the penalty parameter, as a fraction of the maximum value}
 \item{gamma}{Fitting parameter used for tuning backtracking (between 0 and 1)}
\item{step}{Fitting parameter used for initial backtracking step size (between 0 and 1)}
 \item{reset}{Fitting parameter used for taking advantage of local strong convexity in Nesterov momentum (number of iterations before momentum term is reset)}
 \item{cre.sig}{}
 \item{de.sig}{}
 \item{ncores}{Number of computer cores used in calculations}
\item{lambdas}{Optional sequence of lambda values for fitting. We recommend leaving this NULL and letting SGL self-select values}
  \item{verbose}{logical for verbosity of display}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with components 
 \item{best.lambdas}{Best lambda value for each outer fold}
 \item{best.alphas}{Best alpha value for each outer fold}
 \item{pred}{predicted responses (=probabilities for logistic regression)}
 \item{accuracy}{Classification rates}
 \item{ROC}{ROC curve and AUC}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
David Degras and Kourosh Zarringhalam
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

%%\examples{
%%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
