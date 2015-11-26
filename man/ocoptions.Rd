\name{ocoptions}
\alias{ocoptions}
\alias{occolors}
\title{
Options for the opticut package
}
\description{
A convenient way of handling options related to the opticut package.
}
\usage{
ocoptions(...)

occolors(theme)
}
\arguments{
  \item{theme}{
character value, character vector, or a function used to interpolate the colors,
The built-in values are \code{"br"} (blue-red divergent palette, colorblind safe),
\code{"br"} (green-red divergent palette), \code{"bw"} (black and white).
See \code{\link{colorRampPalette}} and the Examples.
}
  \item{\dots}{
arguments in \code{tag = value} form, or a list of tagged values.
The tags must come from the parameters described below.
}
}
\value{
When parameters are set by \code{ocoptions}, their former values are
returned in an invisible named list. Such a list can be passed as an
argument to \code{ocoptions} to restore the parameter values.
Tags are the following:

  \item{collapse}{character value to be used when merging factor levels,
  the default is \code{" "}.}
  \item{cut}{log likelihood ratio value, model/species with lower
  values are excluded from summaries and plots,
  the default is \code{2}.}
  \item{sort}{logical value indicating if species/partitions
  should be meaningfully sorted, the default is \code{TRUE}.}
  \item{theme}{the color theme to be used based on \code{occolors},
  the default is \code{"br"}.}
}
\author{
Peter Solymos <solymos@ualberta.ca>

Hexadecimel values for the built-in palettes are taken from
\url{http://colorbrewer2.org/}.
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

## using palettes
plot(1:100, rep(2, 100), pch = 15,
    ylim = c(0, 15), axes = FALSE, ann = FALSE,
    col = occolors()(100)) # default 'bg'
text(50, 1, "theme = 'br'")
points(1:100, rep(5, 100), pch = 15,
    col=occolors("gr")(100))
text(50, 4, "theme = 'gr'")
points(1:100, rep(8, 100), pch = 15,
    col=occolors("bw")(100))
text(50, 7, "theme = 'bw'")
points(1:100, rep(11, 100), pch = 15,
    col=occolors(terrain.colors)(100))
text(50, 10, "theme = terrain.colors")
points(1:100, rep(14, 100), pch = 15,
    col=occolors(c("orange", "tan", "pink"))(100))
text(50, 13, "theme = c('orange', 'tan', 'pink')")
}
\keyword{ utilities }
