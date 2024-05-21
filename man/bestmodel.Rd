\name{bestmodel}
\alias{bestmodel}
\alias{bestpart}
\alias{getMLE}
\title{
Best model, Partition, and MLE
}
\description{
Generic functions for accessing best model, best partition,
and Maximum Likelihood Estimate from fitted objects.
}
\usage{
bestmodel(object, ...)
bestpart(object, ...)
getMLE(object, ...)
}
\arguments{
  \item{object}{
fitted model object.
}
\item{\dots}{
other arguments passed to the underlying functions.
}
}
\value{
\code{bestmodel} returns the best supported model for further
manipulation (e.g. prediction).

\code{bestpart} returns a matrix with the best supported
partitions for each species (species as columns).

\code{getMLE} returns a named list corresponding to the best supported
model. The list has the following elements:
\code{coef} is the Maximum Likelihood Estimate (MLE),
\code{vcov} is the variance-covariance matrix for the MLE,
\code{dist} is the distribution inherited from input \code{object}.
}
\author{
Peter Solymos <psolymos@gmail.com>
}
\seealso{
\code{\link{opticut}}, \code{\link{multicut}}, \code{\link{uncertainty}}.
}
\keyword{ methods }
