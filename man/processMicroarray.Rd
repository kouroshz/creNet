\name{processMicroarray}
\alias{processMicroarray}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Preprocess MicroArray Data}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
processMicroarray(ents.file, rels.file, celPath, phenoFile, annoFile, ii)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ents.file}{
%%     ~~Describe \code{ents.file} here~~
}
  \item{rels.file}{
%%     ~~Describe \code{rels.file} here~~
}
  \item{celPath}{
%%     ~~Describe \code{celPath} here~~
}
  \item{phenoFile}{
%%     ~~Describe \code{phenoFile} here~~
}
  \item{annoFile}{
%%     ~~Describe \code{annoFile} here~~
}
  \item{ii}{
%%     ~~Describe \code{ii} here~~
}
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ents.file, rels.file, celPath, phenoFile, annoFile, 
    ii) 
{
    L = genOneLevel(ents.file, rels.file)
    ents = L$ents
    rels = L$rels
    ents.mRNA = ents[which(ents$type == "mRNA"), ]
    L = normalizeCEL(ents, rels, celPath, phenoFile, annoFile, 
        ii)
    ents = L$ents
    rels = L$rels
    eset = L$eset
    x.train = L$x.train
    y.train = L$y.train
    x.test = L$x.test
    y.test = L$y.test
    L = list(ents = ents, rels = rels, x.train = x.train, y.train = y.train, 
        x.test = x.test, y.test = y.test)
    L
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
