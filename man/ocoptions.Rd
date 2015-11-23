\name{ocoptions}
\alias{ocoptions}
\title{
Options for the opticut package
}
\description{
A convenient way of handling options related to the opticut package.
}
\usage{
ocoptions(...)
}
\arguments{
  \item{\dots}{
Arguments in \code{tag = value} form, or a list of tagged values.
The tags must come from the parameters described below.
}
}
\value{
When parameters are set by \code{ocoptions}, their former values are 
returned in an invisible named list. Such a list can be passed as an 
argument to \code{ocoptions} to restore the parameter values.
Tags are the following:

  \item{collapse}{Character value to be used when merging factor levels.
  The default is \code{" "}.}
}
\author{
Peter Solymos <solymos@ualberta.ca>
}
\examples{
## resetting pboptions
ocop <- ocoptions(collapse="+")
## check new settings
getOption("ocoptions")
## running again with new settings

## resetting original
ocoptions(opb)
## check reset
getOption("ocoptions")

}
\keyword{ IO }
\keyword{ utilities }
