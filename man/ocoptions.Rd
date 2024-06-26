\name{ocoptions}
\alias{ocoptions}
\title{
Options for the opticut Package
}
\description{
A convenient way of handling options related to the opticut package.
}
\usage{
ocoptions(...)
}
\arguments{
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
  the default is \code{"+"}.}
\item{cut}{log likelihood ratio value, model/species with lower
  values are excluded from summaries and plots,
  the default is \code{2}.}
\item{sort}{logical value indicating if species/partitions
  should be meaningfully sorted, the default is \code{TRUE}.
  It can take numeric value when only species (\code{1})
  or partitions (\code{2}) are to be sorted (\code{1:2} is equivalent to
  \code{TRUE}).}
\item{theme}{the color theme to be used based on \code{\link{occolors}},
  the default is \code{"br"}.}
\item{check_comb}{check the design matrices for complementary partitions
  using \code{\link{checkComb}},
  the default is \code{TRUE}.}
\item{try_error}{if \code{\link{opticut}} and \code{\link{multicut}}
  should \code{\link{try}} to exclude
  species where the models failed (\code{TRUE}), the default
  is to stop when an error is encountered (\code{FALSE}).}
\item{scale}{the scaling factor used to calculate indicator potential (I)
  based on the estimated contrast (x): I = abs(tanh(x*scale)),
  the default is 0.5.}
\item{fix_fitted}{\code{\link{bestpart.multicut}}
  uses \code{\link{lorenz}} which requires nonnegative fitted values,
  however models with identity link can lead to negative expected values.
  When \code{TRUE} the fitted values (x) are adjusted as x' = x + abs(min(x))
  to ensure nonnegativity. The default is \code{FALSE}.}
\item{robust_loglik}{if ill-defined models resulting in perfect fit
  (infinite log likelihood, or \code{NA}, \code{NaN})
  should be allowed. The default \code{TRUE} makes such ill-defined
  log likelihoods a very small real number
  \code{-(.Machine$double.xmax^(1/3))}.
  \code{FALSE} is equivalent to allowing every model to
  safeguard against such cases or not.}
}
\author{
Peter Solymos <psolymos@gmail.com>
}
\examples{
## simple example from Legendre 2013
## Indicator Species: Computation, in
## Encyclopedia of Biodiversity, Volume 4
## https://dx.doi.org/10.1016/B978-0-12-384719-5.00430-5
gr <- as.factor(paste0("X", rep(1:5, each=5)))
spp <- cbind(Species1=rep(c(4,6,5,3,2), each=5),
    Species2=c(rep(c(8,4,6), each=5), 4,4,2, rep(0,7)),
    Species3=rep(c(18,2,0,0,0), each=5))
rownames(spp) <- gr
## must add some noise to avoid perfect fit
spp[6, "Species1"] <- 7
spp[1, "Species3"] <- 17
spp

## current settings
print(unlist(ocoptions())) # these give identical answers
unlist(getOption("ocoptions"))
summary(ocall <- opticut(spp ~ 1, strata=gr, dist="gaussian", comb="all"))

## resetting pboptions and checking new settings
ocop <- ocoptions(collapse="&", sort=FALSE)
unlist(getOption("ocoptions"))
## running again with new settings
summary(ocall <- opticut(spp ~ 1, strata=gr, dist="gaussian", comb="all"))

## resetting original
ocoptions(ocop)
unlist(getOption("ocoptions"))
}
\keyword{ utilities }
