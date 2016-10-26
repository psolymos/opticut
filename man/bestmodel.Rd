\name{bestmodel}
\alias{bestmodel}
\alias{bestpart}
\title{
Best model and partition
}
\description{
Accessing best model and partition from fitted object.
}
\usage{
bestmodel(object, ...)

bestpart(object, ...)
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
}
\author{
Peter Solymos <solymos@ualberta.ca>
}
\seealso{
\code{\link{opticut}}, \code{\link{uncertainty}}.
}
\keyword{ methods }
