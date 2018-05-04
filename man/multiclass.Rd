\name{multiclass}
\alias{multiclass}
\alias{ctable}
\alias{btable}
\alias{untable}
\alias{etable}
\alias{rtable}
\alias{kappacoef}
\alias{test_table}
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
rtable(n, table, method=c("r", "rc"))
kappacoef(predicted, reference)
test_table(table, FUN=kappacoef, n=0,
    type="cohen", method="r", w=NULL, ...)
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
  \item{method}{
method
}
  \item{n}{
n
}
  \item{w}{
w
}
    \item{predicted, reference}{
predicted, reference
}
  \item{FUN}{
FUN
}
  \item{\dots}{
other arguments passed to the \code{FUN}.
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
## contingency table (K x K)
(ct <- ctable(x, y))
## binary contingency tables (2 x 2) by levels
## tp/fp: true/false positives, tn/fn: true/false negatives
(cm <- btable(ct))
## reconstructing the input vectors
untable(ct)
stopifnot(all(ctable(untable(ct)$x, untable(ct)$y) == ctable(x, y)))
## multiclass measures
(x <- multiclass(x, y))
## expected table with same row and col margins
etable(ct, type="cohen")
## single random table based on expected
rtable(1, ct, method="rc")
## compare ct against expected using kappacoef
test_table(ct, FUN=kappacoef, n=0, type="cohen", method="rc")
## compare ct against random distribution
tt <- test_table(ct, FUN=kappacoef, n=99, type="cohen", method="rc")
summary(t(tt))
}
\keyword{ utilities }
