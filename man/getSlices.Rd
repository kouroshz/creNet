\name{getSlices}
\alias{getSlices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Duplicate Covariates According to Regulatory Network and Compute CRE
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getSlices(ents, rels, x, y, cre = c("filter", "weight", "none"), cre.sig = 0.01, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ents}{Data frame of entries typically obtained with \code{\link{processMicroArray}}}
 
  \item{rels}{Data frame of regulator/genes relations typically obtained with \code{\link{processMicroArray}}}
 
  \item{x}{Matrix of predictors (gene expression levels)}

  \item{y}{vector or responses}
 
  \item{cre}{CRE operations: \code{filter} significant regulators, only compute \code{weights} with no filtering, or \code{none} (no CRE)}
  
 \item{cre.sig}{significance level for CRE filtering. Defaults to 0.01}
  
  \item{verbose}{logical flag for verbosity level}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
%%  ~~who you are~~
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
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
