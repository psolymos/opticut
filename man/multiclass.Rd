\name{multiclass}
\alias{multiclass}
\alias{ctable}
\alias{btable}
\alias{untable}
\alias{etable}
\alias{rtable}
\alias{kappa}
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
kappa(predicted, reference)
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
    \item{predicted, reference}{
predicted, reference
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

## here comes etable rtable and more?
}
\keyword{ utilities }
