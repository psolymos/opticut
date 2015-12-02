\name{opticut}
\alias{opticut}
\alias{opticut1}
\alias{plot.opticut}
\alias{wplot}
\alias{wplot.opticut}
\alias{wplot.opticut1}
\alias{print.opticut}
\alias{print.opticut1}
\alias{summary.opticut}
\alias{print.summary.opticut}
\title{
Likelihood based optimal partitioning for indicator species analysis
}
\description{
The functions find the best binary partition for each species based on
model selection. Possibly controlling for modifying/confounding
variables. The general algorithm is described in Kemencei et al. 2014.
}
\usage{
opticut1(Y, X, Z, dist = "gaussian", ...)
opticut(formula, data, strata, dist = "gaussian",
    comb = c("rank", "all"), cl = NULL, ...)

wplot(x, ...)
\method{wplot}{opticut1}(x, cut,
    ylim = c(-1, 1), las=1, ylab = "Model weight * Association",
    xlab = "Partitions", theme, ...)
\method{wplot}{opticut}(x, which = NULL, cut, sort,
    las = 1, ylab = "Model weight * Association",
    xlab = "Partitions", theme, ...)
\method{plot}{opticut}(x, which = NULL, cut, sort,
    las = 1, ylab = "Relative abundance", xlab = "Partitions",
    show_I = TRUE, show_S = TRUE, hr = TRUE,
    theme, mar = c(5, 4, 4, 4) + 0.1, ...)

\method{print}{opticut1}(x, cut, sort, digits, ...)
\method{print}{opticut}(x, digits, ...)
\method{print}{summary.opticut}(x, cut, sort, digits, ...)
\method{summary}{opticut}(object, ...)
}
\arguments{
  \item{formula}{
two sided model formula, response species data (matrix,
or possible a vector for single species case) in the left-hand side,
model terms for modifying effects in the right-hand side
(its structure depending on the underlying functions).
}
  \item{data}{
an optional data frame, list or environment containing the variables
in the model. If not found in data, the variables are taken from
\code{parent.frame()}, typically the environment
from which \code{opticut} is called.
}
  \item{strata}{
vector, unique values define partitions.
}
  \item{dist}{
character or function, a distribution to fit, see Details and Examples.
}
  \item{comb}{
character, how to define the binary partitions.
\code{"rank"} uses \code{\link{rankComb}},
\code{"all"} uses \code{\link{allComb}}.
}
  \item{cl}{
a cluster object, or an integer for multiple cores in parallel computations
(this does nothing on Windows).
}
  \item{Y}{
numeric, vector of observations.
}
  \item{X}{
numeric, design matrix.
}
  \item{Z}{
factor (must have at least 2 unique levels,
this triggers \code{\link{rankComb}}),
or a design matrix (custom matrix or as returned by \code{\link{allComb}}.
}
  \item{x, object}{
object to plot, print, summarize.
}
  \item{cut}{
log likelihood ratio value for showing species.
}
  \item{sort}{
logical, if species should be sorted by indicator values in the output.
}
  \item{show_I}{
logical, if indicator values should be shown.
}
  \item{show_S}{
logical, if number of indicator species should be shown.
}
  \item{hr}{
logical, if horizontal rules should be added.
}
  \item{theme}{
color theme as defined by \code{\link{occolors}}.
}
  \item{mar}{
numeric, graphical parameters for plot margin \code{\link{par}}.
}
  \item{ylab, xlab, las, ylim}{
graphical arguments.
}
  \item{digits}{
numeric, number of significant digits in output.
}
  \item{which}{
numeric, character (can be a vector) defining
a subset of species from the fitted object,
or or \code{NULL} (all species, default).
}
  \item{\dots}{
other arguments passed to the underlying functions, see \code{\link{opticut}}.
}
}
\details{
Currently available distributions:

\code{"gaussian"}: real valued continuous observations, e.g. biomass,
uses \code{\link[stats]{lm}} of the stats package.

\code{"poisson"}: Poisson count data,
uses \code{\link[stats]{glm}} of the stats package.

\code{"binomial"}: presence-absence (detection-nondetection) type data,
uses \code{\link[stats]{glm}} of the stats package.

\code{"negbin"}: overdispersed Negative Binomial count data,
uses \code{\link[MASS]{glm.nb}} of the MASS package.

\code{"beta"}: continuous response in the unit interval, e.g. percent cover,
uses \code{\link[betareg]{betareg}}  of the betareg package.

\code{"zip"}: zero-inflated Poisson counts,
indicative properties are tested as part of the abundance model,
uses \code{\link[pscl]{zeroinfl}} of the pscl package.

\code{"zinb"}: zero-inflated Negative Binomial counts,
indicative properties are tested as part of the abundance model,
uses \code{\link[pscl]{zeroinfl}} of the pscl package.

\code{"zip2"}: zero-inflated Poisson counts,
indicative properties are tested as part of the zero-model,
uses \code{\link[pscl]{zeroinfl}} of the pscl package.

\code{"zinb2"}: zero-inflated Negative Binomial counts,
indicative properties are tested as part of the zero-model,
uses \code{\link[pscl]{zeroinfl}} of the pscl package.

\code{"ordered"}: response measured on ordinal scale,
e.g. ordinal vegetation cover,
uses \code{\link[MASS]{polr}} of the MASS package.

\code{"rsf"}: presence-only data using resource selection
functions (RSF) as explained in  \code{\link[ResourceSelection]{rsf}}
in the ResourceSelection package.
The \code{"rsf"} works only for single species using \code{opticut1}
because 'presence-only' type data cannot be kept in a single
matrix-like object for multiple species.
Intercept only model is accepted for \code{"rsf"}.

\code{"rspf"}: presence-only data using resource selection
probability functions (RSPF)
as explained in \code{\link[ResourceSelection]{rspf}}
in the ResourceSelection package.
The \code{"rspf"} works only for single species using \code{opticut1}
because 'presence-only' type data cannot be kept in a single
matrix-like object for multiple species.
Intercept only model is not accepted for \code{"rspf"}.

Custom distributions can be defined, see Examples.
}
\section{Warning}{
The use of the \code{opticut1} function is generally discouraged.
Use the \code{opticut} function instead with a single species in the left-hand
side of the formula. Some of the internal checks are not guaranteed to
flag issues when the formula-to-model-matrix translation is side-stepped
(this is what is happening when the modifier variables are supplied
as \code{X} argument in \code{opticut1}).
}
\value{
\code{opticut1} returns an object of class opticut1.

\code{opticut} returns an object of class opticut.
}
\references{
Kemencei, Z., Farkas, R., Pall-Gergely, B., Vilisics, F., Nagy, A., Hornung,
E. & Solymos, P. (2014): Microhabitat associations of land snails in
forested dolinas: implications for coarse filter conservation.
Community Ecology 15:180--186.
}
\author{
Peter Solymos <solymos@ualberta.ca>
}
\seealso{
\code{\link{allComb}}, and \code{\link{rankComb}}
for partitioning algorithms.

\code{\link{bestmodel}}, \code{\link{bestpart}}, and \code{\link{uncertainty}}
for manipulating fitted objects.

\code{\link{ocoptions}} on how to set some of the global options
related to the presentation of the results in the package.

\code{\link{optilevel}} for finding the optimal number of factor levels.
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

## all partitions
summary(ocall <- opticut(spp ~ 1, strata=gr, dist="gaussian", comb="all"))

## rank based partitions
summary(ocrank <- opticut(spp ~ 1, strata=gr, dist="gaussian", comb="rank"))

## --- Binomial
## simulated binary data
set.seed(1234)
n <- 200
x0 <- sample(1:4, n, TRUE)
x1 <- ifelse(x0 \%in\% 1:2, 1, 0)
x2 <- rnorm(n, 0.5, 1)
p1 <- plogis(-0.5 + 2*x1 + -0.8*x2)
Y1 <- rbinom(n, 1, p1)
p2 <- plogis(-0.1 + 2*ifelse(x0==4,1,0) + -0.8*x2)
Y2 <- rbinom(n, 1, p2)
p3 <- plogis(-0.1 + -0.8*x2)
Y3 <- rbinom(n, 1, p3)
Y <- cbind(SPP1=Y1, SPP2=Y2, SPP3=Y3)
X <- model.matrix(~x2)

## all partitions, single species
Z <- allComb(x0)
opticut1(Y1, X, Z, dist="binomial")

## rank based partitions, single species
opticut1(Y1, X, as.factor(x0), dist="binomial")

## all partitions, multiple species
m <- opticut(Y ~ x2, strata=x0, dist="poisson", comb="all")
summary(m)
## show all species
summary(m, cut=0)
## model weights
wplot(m)
## associations
plot(m)

## rank based partitions, multiple species
summary(opticut(Y ~ x2, strata=x0, dist="poisson", comb="rank"))

\dontrun{
## --- Zero-inflated Negative Binomial
## dolina analysis
data(dolina)
## stratum as ordinal
dolina$samp$stratum <- as.integer(dolina$samp$stratum)
## filter species to speed up things a bit
Y <- dolina$xtab[,colSums(dolina$xtab > 0) >= 20]
## opticut results
dol <- opticut(Y ~ stratum + lmoist + method, data=dolina$samp,
    strata=dolina$samp$mhab, dist="zinb")
summary(dol)
plot(dol)

## parallel computing comparisons
library(parallel)
cl <- makeCluster(2)
## sequential, all combinations (2^(K-1) - 1)
system.time(opticut(Y ~ stratum + lmoist + method, data=dolina$samp,
    strata=dolina$samp$mhab, dist="zinb", comb="all", cl=NULL))
## sequential, rank based combinations (K - 1)
system.time(opticut(Y ~ stratum + lmoist + method, data=dolina$samp,
    strata=dolina$samp$mhab, dist="zinb", comb="rank", cl=NULL))
## parallel, all combinations (2^(K-1) - 1)
system.time(opticut(Y ~ stratum + lmoist + method, data=dolina$samp,
    strata=dolina$samp$mhab, dist="zinb", comb="all", cl=cl))
## parallel, rank based combinations (K - 1)
system.time(opticut(Y ~ stratum + lmoist + method, data=dolina$samp,
    strata=dolina$samp$mhab, dist="zinb", comb="rank", cl=cl))
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
opticut1(Y[,"amin"], Xdol, Z=dolina$samp$mhab,
    dist=fun, zi_term=dolina$samp$method)
dol2 <- opticut(Y ~ stratum + lmoist + method, data=dolina$samp,
    strata=dolina$samp$mhab, dist=fun, zi_term=dolina$samp$method)
}
}
\keyword{ models }
