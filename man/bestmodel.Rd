\name{bestmodel}
\alias{bestmodel}
\alias{bestmodel.opticut}
\alias{bestmodel.optilevel}
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
\method{bestmodel}{optilevel}(object, ...)

bestpart(object, ...)
\method{bestpart}{opticut}(object, ...)
}
\arguments{
  \item{object}{
fitted model object from \code{\link{opticut}} or
\code{\link{optilevel}}.
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
\code{\link{opticut}} for the user interface.
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

## best partition
head(bestpart(oc))

## best model
mods <- bestmodel(oc)
mods
## explore further
str(predict(mods[[1]]))
confint(mods[[1]])
}
\keyword{ manip }
\keyword{ models }
