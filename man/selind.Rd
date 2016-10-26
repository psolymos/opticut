\name{selind}
\alias{wrsi}
\alias{selind}
\title{
Weighted relative selection index
}
\description{
Calculates weighted relative selection index.
}
\usage{
selind(y, x)
wrsi(y, x)
}
\arguments{
  \item{y}{
matrix of observations for \code{selind},
vector of observations for \code{wrsi}.
}
  \item{x}{
a matrix of proportions (i.e. the values 0 and 1 should have
consistent meaning across the columns, often through a unit sum constraint).
}
}
\value{
\code{wrsi} returns a data frame (class 'wrsi') with the following columns:
\describe{
\item{\code{WRSI}}{weighted relative selection index, range (0- Inf).}
\item{\code{zWRSI}}{log of \code{WRSI} (z-transformed), range (-Inf, Inf).}
\item{\code{rWRSI}}{inverse Fisher z-transformed \code{zWRSI}, range (-1, 1).}
\item{\code{Pused} and \code{Pavail}}{total proportion of used (\code{y > 0})
and available of each feature (column) in \code{x}.}
\item{\code{Pw}}{weighted proportions from \code{y}.}
\item{\code{u} and \code{a}}{used and available totals
for each feature (column) in \code{x}.}
}

\code{selind} returns a data frame (class 'selind') with
one column for each species, and one row for each feature (column) in \code{x}.
Cell values are inverse Fisher z-transformed (\code{zWRSI}) indices.
}
\references{
Authors
Lele, S. R., Merrill, E. H., Keim, J., Boyce, M. S., 2013.
Selection, use, choice and occupancy:
clarifying concepts in resource selection studies.
Journal of Animal Ecology 82:1183-1191.

}
\author{
Peter Solymos <solymos@ualberta.ca>
}
\seealso{
\code{\link{optilevels}} for finding optimal number of factor levels.
}
\examples{
## --- habitat composition matrix
set.seed(1234)
n <- 1000 # sample size
k <- 5 # habitat classes
s <- runif(n, 1, 5)
p <- plogis(rnorm(n*k, 0, rep(s, k)))
p <- p*t(replicate(n, sample(c(10,4,2,1,1))))
x <- p / rowSums(p)
summary(x)
summary(rowSums(x))

## --- observations
## expected abundance in each habitat class
lam <- c(0.8, 0.6, 0.5, 0.4, 0.1)*1
## sample x habitat level abundances
yy <- t(sapply(seq_len(n), function(i) {
    ## intercept and modifier combined
    rpois(k, (x[i,]*lam))
    }))
## total: sum over habitat classes
## this is what we observe
y <- rowSums(yy)
colSums(yy)
table(y)

## --- wrsi calculations
(w <- wrsi(y, x))
op <- par(mfrow=c(1,2))
## habitat level observations are unknown
plot(lam, colSums(yy) / sum(yy), type="b")
## this is approximated by the wrsi
plot(lam, w$rWRSI, type="b")
abline(h=0, lty=2)
par(op)

## --- selind calculations for multiple species
y2 <- cbind(Spp1=y, Spp2=rev(y), Spp3=sample(y))
(w2 <- selind(y2, x))
heatmap(as.matrix(w2), scale="none", col=occolors()(25))
}
\keyword{ models }
