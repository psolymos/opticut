\name{allComb}
\alias{allComb}
\alias{kComb}
\alias{checkComb}
\title{
Finding All Possible Binary Partitions
}
\description{
These functions are used to find all possible binary partitions.
Finding all combinations require
a classification vector with K > 1 strata.
}
\usage{
allComb(x, collapse)
kComb(k)
checkComb(x)
}
\arguments{
  \item{x}{
a vector for \code{allComb} (can be of any type but treated as factor, must have at least 2 unique values); and a numeric matrix for \code{checkComb}.
}
  \item{collapse}{
character, what to paste between levels.
Defaults to \code{getOption("ocoptions")$collapse}.
}
  \item{k}{
numeric, number of levels (strata) in a given classification (K > 1).
}
}
\value{
\code{kComb} returns a contrast matrix corresponding to
all possible binary partitions of the factor with K levels.
Complements are not counted twice, i.e.
(0,0,1,1) is equivalent to (1,1,0,0).
The number of such possible combinations is M = 2^(K - 1) - 1.

\code{allComb} takes a classification vector with at least 2 levels
and returns a model matrix with binary partitions.

\code{checkComb} checks if combinations are unique and non-complementary
(misfits are returned as attributes). Returns a logical value.
}
\author{
Peter Solymos <psolymos@gmail.com>
}
\seealso{
\code{\link{opticut}} for the user interface.

\code{\link{rankComb}} and \code{\link{lorenz}}
for alternative partitioning algorithms.
}
\examples{
kComb(k = 2)
kComb(k = 3)
kComb(k = 4)

## finding all combinations
(f <- rep(LETTERS[1:4], each=2))
(mc <- allComb(f, collapse = "_"))
## checking for complementary entries
checkComb(mc) # TRUE
## adding complementary entries to the matrix
mc2 <- cbind(z = 1 - mc[,1], mc[,c(1:ncol(mc), 1)])
colnames(mc2) <- 1:ncol(mc2)
mc2
checkComb(mc2) # FALSE
}
\keyword{ manip }
\keyword{ misc }
