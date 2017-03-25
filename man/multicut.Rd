\name{multicut}
\alias{multicut}
\alias{multicut.default}
\alias{multicut.formula}
\alias{multicut1}
\alias{getMLE.multicut}
\alias{bestpart.multicut}
\alias{bestmodel.multicut}
\alias{fitted.multicut}
\alias{predict.multicut}
\alias{plot.multicut}
\alias{strata.multicut}
\alias{print.multicut}
\alias{print.multicut1}
\alias{subset.multicut}
\alias{summary.multicut}
\alias{print.summary.multicut}
\alias{as.data.frame.multicut}
\alias{as.data.frame.summary.multicut}
\alias{getMLE.multicut}
\title{
Multi-level Response Model
}
\description{
The functions fits the multi-level response model for each species,
possibly controlling for modifying/confounding variables.
}
\usage{
multicut1(Y, X, Z, dist = "gaussian", sset=NULL, ...)

multicut(...)
\method{multicut}{default}(Y, X, strata, dist = "gaussian",
    sset=NULL, cl = NULL, ...)
\method{multicut}{formula}(formula, data, strata, dist = "gaussian",
    sset=NULL, cl = NULL, ...)

\method{bestmodel}{multicut}(object, which = NULL, ...)
\method{bestpart}{multicut}(object, ...)
\method{strata}{multicut}(object, ...)
\method{getMLE}{multicut}(object, which, vcov=FALSE, ...)
\method{subset}{multicut}(x, subset=NULL, ...)
\method{fitted}{multicut}(object, ...)
\method{predict}{multicut}(object, gnew=NULL, xnew=NULL, ...)

\method{plot}{multicut}(x, which = NULL, cut, sort,
    las, ylab = "Relative abundance", xlab = "Strata",
    show_I = TRUE, show_S = TRUE, hr = TRUE, tick = TRUE,
    theme, mar = c(5, 4, 4, 4) + 0.1, bty = "o",
    lower = 0, upper = 1, pos = 0, horizontal=TRUE, ...)

\method{print}{multicut1}(x, digits, ...)
\method{print}{multicut}(x, digits, ...)
\method{print}{summary.multicut}(x, cut, sort, digits, ...)
\method{summary}{multicut}(object, ...)

\method{as.data.frame}{multicut}(x,
    row.names = NULL, optional = FALSE, cut, sort, ...)
\method{as.data.frame}{summary.multicut}(x,
    row.names = NULL, optional = FALSE, cut, sort, ...)
}
\arguments{
  \item{formula}{
two sided model formula, response species data (matrix,
or possible a vector for single species case) in the left-hand side,
model terms for modifying effects in the right-hand side
(its structure depending on the underlying functions).
For example, in the most basic Gaussian case it can be
\code{y ~ 1} (no modifying variables) or \code{y ~ x}
(with modifying variables).
Centering the modifying terms (or choosing the origin wisely)
is generally recommended (especially for Gaussian distribution
where linear predictors are additive on the response scale)
because the relative abundance contrasts are estimated at the origin (0).
}
  \item{data}{
an optional data frame, list or environment containing the variables
in the model. If not found in data, the variables are taken from
\code{parent.frame()}, typically the environment
from which \code{multicut} is called.
}
  \item{strata}{
a factor, unique values define partitions
(must have at least 2 unique levels, empty levels are dropped).
}
  \item{dist}{
character or function, a distribution to fit.
If character, it can follow one of these patterns: \code{"family"},
or \code{"family:link"} when appropriate (there is a \code{link}
argument in the underlying function, or the link
can be specified via the \code{family} argument).
See Details on \code{\link{opticut}} page and Examples.
}
  \item{sset}{
an optional vector specifying a subset of observations (rows)
to be used in the fitting process. \code{NULL} means no subset taken.
}
  \item{cl}{
a cluster object, or an integer for multiple cores in parallel computations
(integer value for forking is ignored on Windows).
}
  \item{Y}{
numeric vector of observations for \code{multicut1},
vector or community matrix for \code{multicut.default}.
}
  \item{X}{
numeric, design matrix. Can be missing, in which case an intercept-only model
is assumed.
}
  \item{Z}{
factor (must have at least 2 unique levels).
}
  \item{x, object}{
object to plot, print, summarize.
}
  \item{cut}{
log likelihood ratio value to be used as a cut-off for showing species
whose log likelihood ratio is not less than the cut-off.
}
  \item{sort}{
logical value indicating if species/partitions
should be meaningfully sorted, the default is \code{TRUE}.
It can take numeric value when only species (\code{1})
or partitions (\code{2}) are to be sorted (\code{1:2} is equivalent to
\code{TRUE}).
}
  \item{show_I}{
logical, if indicator potential (I) should be shown.
}
  \item{show_S}{
logical, if number of indicator species should be shown.
}
  \item{hr, tick}{
logical, if horizontal rules (\code{hr})
and ticks to the axis legends (\code{tick}) should be added.
Default is \code{TRUE} for both.
}
  \item{theme}{
color theme as defined by \code{\link{occolors}}.
}
  \item{mar}{
numeric, graphical parameters for plot margin \code{\link{par}}.
}
  \item{ylab, xlab, las}{
graphical arguments, see \code{\link{plot}}.
By default, \code{las} is 1 when \code{horizontal = TRUE} and
2 when \code{horizontal = FALSE}.
}
  \item{bty}{
Character, determines the type of box which is drawn around plots,
see \code{\link{par}}.
}
  \item{lower, upper}{
numeric (between 0 and 1), \code{lower} is the minimum and
\code{upper} is the maximum height for rectangles drawn in the plot.
Both need to be in [0, 1] and \code{higher} cannot be smaller than \code{lower}.
}
  \item{pos}{
numeric, position of rectangles in the plot relative to the baseline.
Value must be in the [-1, 1] range (below vs. above baseline).
}
  \item{horizontal}{
logical, plot orientation: species as rows (\code{TRUE})
or as columns (\code{FALSE}).
}
  \item{digits}{
numeric, number of significant digits in output.
}
  \item{which}{
numeric or character (can be a vector) defining
a subset of species from the fitted object,
or \code{NULL} (all species, default).
}
\item{row.names}{
\code{NULL} or a character vector giving the row names for the data frame.
Missing values are not allowed. See \code{\link{as.data.frame}}.
}
\item{optional}{
logical. If \code{TRUE}, setting row names and converting column names
(to syntactic names: see \code{\link{make.names}}) is optional.
See \code{\link{as.data.frame}}.
}
  \item{subset}{
logical, numeric, or character index indicating species to keep,
missing values are not accepted. The default \code{NULL} returns
the original object without subsetting.
}
  \item{vcov}{
logical, if variance-covariance matrix is to be returned.
}
  \item{gnew, xnew}{
new values for \code{strata} and modifiers (right-hand-side of formula)
to predict for, or \code{NULL}.
}
  \item{\dots}{
other arguments passed to the underlying functions.
}
}
\details{
FIXME

add Lc approach and how I is calculated
}
\section{Warning}{
The use of the \code{multicut1} function is generally discouraged:
some of the internal checks are not guaranteed to
flag issues when the formula-to-model-matrix translation is side-stepped
(this is what is happening when the modifier variables are supplied
as \code{X} argument in \code{multicut1}).
Use the \code{multicut} function with a single species instead.
}
\value{
\code{multicut1} returns an object of class 'multicut1'.

\code{multicut} returns an object of class 'multicut', that is a list
with the following components:
\describe{
\item{\code{"call"}}{the function call.}
\item{\code{"species"}}{a list of species specific multicut1 objects.}
\item{\code{"X"}}{modifying variables as model matrix.}
\item{\code{"Y"}}{response, single species vector or matrix.}
\item{\code{"strata"}}{defines the stratification.}
\item{\code{"nobs"}}{sample size.}
\item{\code{"sset"}}{subset, if specified.}
\item{\code{"dist"}}{distribution.}
\item{\code{"failed"}}{IDs for failed species models dropped from results list.}
}

The \code{strata} method extracts the \code{strata} argument
as factor.

The \code{print} and \code{summary} methods are called for their side effects
showing expected values, and log likelihood ratio (logLR).

The \code{subset} method subsets the species in the multicut object.

The \code{plot} method presents the estimates by species and strata.

\code{\link{bestpart}} returns a matrix with the best supported
partitions for each species (samples and rows, species as columns).
Binary partitions are based on Lorenz-tangent based optimal threshold
(see \code{\link{lorenz}}).
\code{\link{lorenz}} requires nonnegative fitted values
which is not guaranteed under \code{dist = "gaussian"} with identity
link, see \code{fix_fitted} \code{\link{ocoptions}} setting
for how to resolve this (choosing a different link function,
distribution, or centering modified variables is advised).

\code{\link{bestmodel}} returns the best supported model for further
manipulation (e.g. prediction). Note: custom distribution
functions are designed to return only point estimates,
thus the best model cannot be returned. In this case,
use the best partition returned by \code{bestpart} to refit the model.
\code{getMLE} returns a named list corresponding to the best supported
model. The list has the following elements:
\code{coef} is the Maximum Likelihood Estimate (MLE),
\code{vcov} is the variance-covariance matrix for the MLE or \code{NULL},
\code{dist} is the distribution inherited from input \code{object}.

\code{\link{fitted}} returns expected values on the predictor scale
for the observations as a matrix (number of observations by number of species).
\code{\link{predict}} returns \code{fitted} values when both \code{gnew}
and \code{xnew} are \code{NULL}, or corresponding point predictions
(expected values) on the predictor scale.

The coercion methods \code{\link{as.data.frame}} return a data frame.
}
\author{
Peter Solymos <solymos@ualberta.ca>
}
\seealso{
\code{\link{lorenz}} Examples for how multi-level partitions are binarized
using the Lorenz-tangent approach.

\code{\link{opticut}} for optimal binary response model,
\code{\link{optilevels}} for finding the optimal number of factor levels.

\code{\link{beta2i}} for indicator potential (I) calculations in summaries.

\code{\link{bestmodel}}, \code{\link{bestpart}}, and \code{\link{uncertainty}}
for manipulating fitted objects.

\code{\link{ocoptions}} on how to set some of the global options
related to the presentation of the results in the package
and how errors encountered during model fitting are handled.
}
\examples{
## --- Gaussian
## simple example from Legendre 2013
## Indicator Species: Computation, in
## Encyclopedia of Biodiversity, Volume 4
## http://dx.doi.org/10.1016/B978-0-12-384719-5.00430-5
gr <- as.factor(paste0("X", rep(1:5, each=5)))
spp <- cbind(Species1=rep(c(4,6,5,3,2), each=5),
    Species2=c(rep(c(8,4,6), each=5), 4,4,2, rep(0,7)),
    Species3=rep(c(18,2,0,0,0), each=5))
rownames(spp) <- gr
spp

## negative expected values are not good
oco <- ocoptions(fix_fitted=TRUE)
summary(ocall <- multicut(spp ~ 1, strata=gr, dist="gaussian"))
summary(multicut(spp, strata=gr, dist="gaussian")) # alternative
ocoptions(oco) # reset options

## --- Binomial
## simulated binary data
set.seed(1234)
n <- 200
x0 <- sample(1:4, n, TRUE)
x1 <- ifelse(x0 <= 2, 1, 0)
x2 <- rnorm(n, 0.5, 1)
p1 <- plogis(-0.5 + 2*x1 + -0.8*x2)
Y1 <- rbinom(n, 1, p1)
p2 <- plogis(-0.1 + 2*ifelse(x0==4,1,0) + -0.8*x2)
Y2 <- rbinom(n, 1, p2)
p3 <- plogis(-0.1 + -0.8*x2)
Y3 <- rbinom(n, 1, p3)
Y <- cbind(SPP1=Y1, SPP2=Y2, SPP3=Y3)
X <- model.matrix(~x2)

multicut1(Y1, X, as.factor(x0), dist="binomial")

summary(m1 <- multicut(Y ~ x2, strata=x0, dist="poisson"))
plot(m1)

## subset results
summary(subset(m1, 1:2))

## best partition
head(bestpart(m1))

## best model
mods <- bestmodel(m1)
mods
## explore further
confint(mods[[1]])

## MLE and variance-covariance matrix (species 1)
getMLE(m1, which = 1, vcov=TRUE)

## fitted values
head(fitted(m1))
## prediction for new data
head(predict(m1, gnew=x0, xnew=data.frame(x2=x2)))

\dontrun{
## --- Zero-inflated Negative Binomial
## dolina example
data(dolina)
## stratum as ordinal
dolina$samp$stratum <- as.integer(dolina$samp$stratum)
## filter species to speed up things a bit
Y <- dolina$xtab[,colSums(dolina$xtab > 0) >= 20]
## opticut results, note the cloglog link function
dol <- multicut(Y ~ stratum + lmoist + method, data=dolina$samp,
    strata=dolina$samp$mhab, dist="zinb:cloglog")
summary(dol)
## vertical plot orientation
plot(dol, horizontal=FALSE, pos=1, upper=0.8)

## parallel
library(parallel)
cl <- makeCluster(2)
multicut(Y ~ stratum + lmoist + method, data=dolina$samp,
    strata=dolina$samp$mhab, dist="zip",cl=cl)
stopCluster(cl)

## --- Customizing distributions
## we may want to expand the Zero-inflation component in a ZIP model
## see how the return value needs to be structured
fun <- function(Y, X, linkinv, zi_term, ...) {
    X <- as.matrix(X)
    mod <- pscl::zeroinfl(Y ~ X-1 | zi_term, dist = "poisson", ...)
    list(coef=coef(mod),
        logLik=logLik(mod),
        linkinv=mod$linkinv)
}
Xdol <- model.matrix(~ stratum + lmoist + method, data=dolina$samp)
## this fits the null model (i.e. no partitions added)
fun(Y[,"amin"], Xdol, zi_term=dolina$samp$method)
## now we can use dist=fun
multicut1(Y[,"amin"], Xdol, Z=dolina$samp$mhab,
    dist=fun, zi_term=dolina$samp$method)
dol2 <- multicut(Y ~ stratum + lmoist + method, data=dolina$samp,
    strata=dolina$samp$mhab, dist=fun, zi_term=dolina$samp$method)
summary(dol2)
}
}
\keyword{ models }
\keyword{ utilities }
