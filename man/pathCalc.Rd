\name{pathCalc}
\alias{pathCalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select Path of Regularization Parameters 
}
\description{
}
\usage{
pathCalc(data, index, weights, alpha = 0.95, min.frac = 0.05, nlam = 100, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
list with components $x$ and $y$}
  \item{index}{
group membership indicator}
  \item{weights}{
group weights for regularization penalty 
}
  \item{alpha}{
tradeoff between lasso penalty and group lasso penalty. \code{alpha=1} corresponds to pure lasso, \code{alpha=0} to pure group lasso. 
}
  \item{min.frac}{
minimum regularization parameter lambda as a fraction of the largest
}
  \item{nlam}{
number of lambda values along the path}
  \item{type}{
type of regression model: \code{linear} or \code{logit}}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of \code{nlam} lambda values
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
David Degras
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
