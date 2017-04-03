\name{multiclass}
\alias{multiclass}
\alias{ctable}
\alias{btable}
\alias{untable}
\alias{etable}
\alias{rtable}
\title{
Multiclass Measures and Support Functions
}
\description{
Multiclass measures of predictive performance and supporting functions.
}
\usage{
ctable(x, y)
btable(table)
untable(table)

multiclass(x, y=NULL, beta=1)

etable(table, type="cohen", w=NULL)
rtable(n, table, type=c("r", "rc"))
}
\arguments{
  \item{x}{
x
}
  \item{y}{
y
}
  \item{table}{
table
}
  \item{beta}{
beta
}
  \item{type}{
type
}
  \item{n}{
n
}
  \item{w}{
w
}
}
\details{
DETAILS
}
\value{
VALUE
}
\author{
Peter Solymos <solymos@ualberta.ca>
}
\references{
REFS
}
\seealso{
SEEALSO
}
\examples{
N <- 100
K <- 5
x <- sample(LETTERS[1:K], N, replace=TRUE, prob=sqrt(2^(0:(K-1))))
y <- x
for (i in 1:length(x))
    if (runif(1) > 0.1) # error rate
        y[i] <- sample(LETTERS[1:K], 1)
(ct <- ctable(x, y))
(cm <- btable(ct))
untable(ct)
stopifnot(all(ctable(untable(ct)$x, untable(ct)$y) == ctable(x, y)))
(x <- multiclass(x, y))
}
\keyword{ utilities }
