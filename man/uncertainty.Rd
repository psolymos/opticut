\name{uncertainty}
\alias{uncertainty}
\alias{uncertainty.opticut}
\alias{print.uncertainty1}
\alias{print.uncertainty}
\alias{print.summary.uncertainty}
\alias{summary.uncertainty}
\alias{bestpart.uncertainty}
\alias{bestpart.uncertainty1}
\alias{check_strata}
\alias{as.matrix.uncertainty}
\alias{as.matrix.summary.uncertainty}
\alias{as.data.frame.uncertainty}
\alias{as.data.frame.uncertainty}
\title{
Quantifying uncertainty for fitted objects
}
\description{
Quantifying uncertainty for fitted objects.
}
\usage{
uncertainty(object, ...)
\method{uncertainty}{opticut}(object,
    which = NULL, type = c("asymp", "boot", "multi"),
    B = 99, cl = NULL, ...)

check_strata(x, mat)

\method{bestpart}{uncertainty}(object, ...)
\method{bestpart}{uncertainty1}(object, ...)

\method{print}{uncertainty1}(x, ...)
\method{print}{uncertainty}(x, ...)
\method{print}{summary.uncertainty}(x, sort, digits, ...)
\method{summary}{uncertainty}(object, level = 0.95, ...)

\method{as.matrix}{uncertainty}(x, sort, ...)
\method{as.matrix}{summary.uncertainty}(x, sort, ...)
\method{as.data.frame}{uncertainty}(x, sort, ...)
\method{as.data.frame}{summary.uncertainty}(x, sort, ...)
}
\arguments{
  \item{object}{
fitted model object from \code{\link{opticut}} (or an output from
\code{uncertainty} for the \code{summary} method).
}
  \item{which}{
numeric or character (can be a vector) defining
a subset of species from the fitted object,
or or \code{NULL} (all species, default).
}
  \item{type}{
character, describing the type of uncertainty calculation for
indicator value \code{I}, and expected values \code{mu0}, and \code{mu1}.
\code{"asymp"}: asymptotic distribution is based on best partition
found for the input object.
\code{"boot"}: non-parametric bootstrap distribution
based on best partition found for the input object.
\code{"multi"}: non-parametric bootstrap distribution
based on best partition found for the bootstrap data (i.e.
the model ranking is re-evaluated each time).
\code{"multi"} works only if \code{comb = "rank"} in the
\code{\link{opticut}} call.
}
  \item{B}{
numeric, number of iterations. For \code{type = "boot"} and
\code{type = "multi"} it can be a user-supplied matrix with indices for
resampling with dimensions length of observations times B.
}
  \item{cl}{
a cluster object, or an integer for multiple cores in parallel computations
(this does nothing on Windows).
}
  \item{x}{
an object to be printed.
}
  \item{level}{
the confidence level required.
}
  \item{digits}{
numeric, number of significant digits in output.
}
  \item{mat}{
a matrix with resampling indices (rows as samples, columns as iterations).
}
\item{\dots}{
other arguments passed to the underlying functions.
}
}
\value{
\code{uncertainty} returns an object of class 'uncertainty'.
The \code{uncertainty} element of the object is a list with species specific
output as elements (object class 'uncertainty1').
Each 'uncertainty1' output is a data frame with 4 columns:
\code{best} partition, indicator value \code{I},
and expected values \code{mu0}, and \code{mu1}.

\code{check_strata} returns a logical vector checking if
all original strata from the input object are represented
by resampling indices. Number of strata are attached as attributes
for further diagnostics.

The summary method prints the name of the best supported split,
selection frequency (R, reliability), indicator values (I, based on
the distribution of values within the best supported split with highest
reliability) and conficence interval for I (based on \code{level}).
}
\section{Warning}{
Resampling methods can lead to complete exclusion of certain strata
when sample size is small. Try revising the stratification of the
input object, or provide custom resampling indices via the \code{B}
argument using stratified (block) bootstrap, jackknife (leave-one-out),
or similar techniques. Sometimes finding a suitable random seed
via \code{\link{set.seed}} can resolve the issue.
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
x3 <- ifelse(x0 \%in\% 2:4, 1, 0)
lam1 <- exp(0.5 + 1*x1 + -0.2*x2)
Y1 <- rpois(n, lam1)
lam2 <- exp(1 + 0.5*x3)
Y2 <- rpois(n, lam2)
Y3 <- rpois(n, exp(0))
Y <- cbind(Spp1=Y1, Spp2=Y2, Spp3=Y3)

oc <- opticut(Y ~ x2, strata=x0, dist="poisson", comb="rank")

## asymptotic confidence intervals
uc1 <- uncertainty(oc, type="asymp", B=999)
summary(uc1)
## bootstrap-based confidence intervals
uc2 <- uncertainty(oc, type="boot", B=19)
summary(uc2)

## use user-supplied indices
## multi-model bootstrap based uncertainties
B <- replicate(25, sample.int(n, replace=TRUE))
check_strata(oc, B) # check representation
uc3 <- uncertainty(oc, type="multi", B=B)
summary(uc3)

## best partitions:
## selection frequencies for strata and species
bestpart(uc3)
heatmap(bestpart(uc3), scale="none", col=occolors()(25))

## individual species results
uc3$uncertainty
bestpart(uc3$uncertainty[[1]])

## block bootstrap
block_fun <- function()
    unlist(lapply(unique(x0), function(z) if (sum(x0==z) < 2)
        which(x0==z) else sample(which(x0==z), sum(x0==z), replace=TRUE)))
B <- replicate(25, block_fun())
check_strata(oc, B) # check representation
summary(uncertainty(oc, type="multi", B=B))

## jackknife
B <- sapply(1:n, function(i) which((1:n) != i))
check_strata(oc, B) # check representation
summary(uncertainty(oc, type="multi", B=B))

\dontrun{
## dolina example
data(dolina)
## stratum as ordinal
dolina$samp$stratum <- as.integer(dolina$samp$stratum)
## filter species to speed up things a bit
Y <- ifelse(dolina$xtab[,colSums(dolina$xtab > 0) >= 20] > 0, 1, 0)
## opticut results, note the cloglog link function
dol <- opticut(Y ~ stratum + lmoist + method, data=dolina$samp,
    strata=dolina$samp$mhab, dist="binomial:cloglog")

## parallel computing for uncertainty
library(parallel)
cl <- makeCluster(2)
ucdol <- uncertainty(dol, type="multi", B=25, cl=cl)
stopCluster(cl)

bestpart(ucdol)
heatmap(t(bestpart(ucdol)), scale="none", col=occolors()(25),
    distfun=function(x) dist(x, "manhattan"))
}
}
\keyword{ manip }
\keyword{ models }
