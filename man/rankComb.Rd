\name{rankComb}
\alias{rankComb}
\alias{oComb}
\title{
Ranking based binary partitions
}
\description{
Blindly fitting a model to all possible partitions is wasteful
use of resources. Instead, one can rank the K partitions
based on expected response values
to explore only K-1 binary partitions along the gradient
defined by the ranks of the expected values.
}
\usage{
oComb(x, collapse = getOption("ocoptions")$collapse)
rankComb(Y, X, Z, dist = "gaussian", 
    collapse = getOption("ocoptions")$collapse, ...)
}
\arguments{
  \item{Y}{
numeric, vector of observations.
}
  \item{X}{
numeric, design matrix.
}
  \item{Z}{
factor, must have at least 2 unique levels.
}
  \item{dist}{
character, distribution argument passed to underlying functions,
see listed on the help page of \code{\link{opticut}}.
}
  \item{x}{
a vector for \code{allComb} (can be of any type but treated as factor, must have at least 2 unique values); and a numeric matrix for \code{checkComb}.
}
  \item{collapse}{
character, what to paste between levels.
}
  \item{\dots}{
other arguments passed to the underlying functions, see \code{\link{opticut}}.
}
}
\value{
\code{oComb} return the 'contrast' matrix based on the rank vector as input.
Ranked from lowest to highest expected value among the partitions.

The function \code{rankComb} fits the model with multiple (K > 2) factor levels
to find out the ranking, and returns a binary classification matrix.
}
\author{
Peter Solymos <solymos@ualberta.ca>
}
\seealso{
\code{\link{allComb}} for alternative partitioning algorithm.

\code{\link{opticut}} for the user interface.
}
\examples{
## simulate some data
set.seed(1234)
n <- 200
x0 <- sample(1:4, n, TRUE)
x1 <- ifelse(x0 \%in\% 1:2, 1, 0)
x2 <- rnorm(n, 0.5, 1)
lam <- exp(0.5 + 0.5*x1 + -0.2*x2)
Y <- rpois(n, lam)

## binary partitions
head(rc <- rankComb(Y, model.matrix(~x2), as.factor(x0), dist="poisson"))
attr(rc, "est") # expected values in factor levels

oComb(1:4)
oComb(attr(rc, "est"))
}
\keyword{ manip }
\keyword{ misc }
