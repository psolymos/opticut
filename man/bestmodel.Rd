\name{bestmodel}
\alias{bestmodel}
\alias{bestmodel.opticut}
\alias{bestpart}
\alias{bestpart.opticut}
\title{
Accessing best model and partition from fitted object
}
\description{
Accessing best model and partition from fitted object.
}
\usage{
bestmodel(object, ...)
\method{bestmodel}{opticut}(object, which = NULL, ...)
bestpart(object, ...)
\method{bestpart}{opticut}(object, ...)
}
\arguments{
  \item{object}{
fitted model object from \link{\code{opticut}}.
}
  \item{which}{
numeric, character (can be a vector) defining
a subset of species from the fitted object,
or or \code{NULL} (all species, default).
}
  \item{\dots}{
other arguments passed to the underlying functions.
}
}
\value{
\code{bestmodel} returns the best supported mode for further
manipulation (e.g. prediction).

\code{bestpart} returns a matrix with the best supported
partitions for each species (species as columns).
}
\author{
Peter Solymos <solymos@ualberta.ca>
}
\seealso{
\link{\code{opticut}} for the user interface.
}
\examples{
set.seed(2345)
n <- 50
x0 <- sample(1:4, n, TRUE)
x1 <- ifelse(x0 \%in\% 1:2, 1, 0)
x2 <- rnorm(n, 0.5, 1)
lam <- exp(0.5 + 1*x1 + -0.2*x2)
Y <- rpois(n, lam)
oc <- opticut(Y ~ x2, strata=x0, dist="poisson", comb="rank")

head(bestpart(oc))

bestmodel(oc)
}
\keyword{ manip }
\keyword{ models }
