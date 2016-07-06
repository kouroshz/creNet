\name{standardize}
\alias{standardize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standardize A Matrix}
\description{
This function centers and scales the columns of a matrix. It essentially produces the same result as \code{\link[base]{scale}} 
}
\usage{
standardize(x)
}

\arguments{
  \item{x}{a matrix}
}

\value{
A list with components
\item{x}{the standardized matrix \code{x}}
\item{center}{column means of \code{x}}
\item{scale}{norms of the (centered) columns of \code{x}}
}



\examples{
x <- matrix(runif(40),8,5)
x.std <- standardize(x)
colMeans(x.std$x)
colSums(x.std$x^2)
}

