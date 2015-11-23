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
## simple example from Legendre 2013
## Indicator Species: Computation, in
## Encyclopedia of Biodiversity, Volume 4
## http://dx.doi.org/10.1016/B978-0-12-384719-5.00430-5
gr <- as.factor(paste0("X", rep(1:5, each=5)))
spp <- cbind(Species1=rep(c(4,6,5,3,2), each=5),
    Species2=c(rep(c(8,4,6), each=5), 4,4,2, rep(0,7)),
    Species3=rep(c(18,2,0,0,0), each=5))
rownames(spp) <- gr

## current settings
print(unlist(ocoptions())) # these give identical answers
unlist(getOption("ocoptions"))
summary(ocall <- opticut(spp ~ 1, strata=gr, dist="gaussian", comb="all"))


## resetting pboptions and checking new settings
ocop <- ocoptions(collapse="+", sort=FALSE)
unlist(getOption("ocoptions"))
## running again with new settings
summary(ocall <- opticut(spp ~ 1, strata=gr, dist="gaussian", comb="all"))

## resetting original
ocoptions(ocop)
unlist(getOption("ocoptions"))
}
\keyword{ utilities }
