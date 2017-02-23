\name{opticut}
\alias{opticut}
\alias{opticut.default}
\alias{opticut.formula}
\alias{opticut1}
\alias{fitted.opticut}
\alias{predict.opticut}
\alias{plot.opticut}
\alias{wplot}
\alias{fix_levels}
\alias{strata}
\alias{strata.opticut}
\alias{wplot.opticut}
\alias{wplot.opticut1}
\alias{print.opticut}
\alias{print.opticut1}
\alias{summary.opticut}
\alias{print.summary.opticut}
\alias{bestmodel.opticut}
\alias{bestpart.opticut}
\alias{as.data.frame.opticut}
\alias{as.data.frame.summary.opticut}
\alias{getMLE.opticut}
\title{
Likelihood based optimal partitioning for indicator species analysis
}
\description{
The functions find the best binary partition for each species based on
model selection. Possibly controlling for modifying/confounding
variables. The general algorithm is described in Kemencei et al. 2014.
}
\usage{
opticut1(Y, X, Z, dist = "gaussian", sset=NULL, ...)

opticut(...)
\method{opticut}{default}(Y, X, strata, dist = "gaussian",
    comb = c("rank", "all"), sset=NULL, cl = NULL, ...)
\method{opticut}{formula}(formula, data, strata, dist = "gaussian",
    comb = c("rank", "all"), sset=NULL, cl = NULL, ...)

fix_levels(x, sep = "_")
strata(object, ...)
\method{strata}{opticut}(object, ...)

\method{bestmodel}{opticut}(object, which = NULL, ...)
\method{bestpart}{opticut}(object, pos_only = FALSE, ...)
\method{getMLE}{opticut}(object, which, ...)
\method{fitted}{opticut}(object, ...)
\method{predict}{opticut}(object, gnew=NULL, xnew=NULL, ...)

wplot(x, ...)
\method{wplot}{opticut1}(x, cut, ylim = c(-1, 1),
    las=1, ylab = "Model weight * Association", xlab = "Partitions",
    theme, mar = c(5, 4, 4, 4) + 0.1, bty = "o", ...)
\method{wplot}{opticut}(x, which = NULL, cut, sort,
    las = 1, ylab = "Model weight * Association", xlab = "Partitions",
    theme, mar = c(5, 4, 4, 4) + 0.1, bty = "o", ...)
\method{plot}{opticut}(x, which = NULL, cut, sort,
    las, ylab = "Relative abundance", xlab = "Strata",
    show_I = TRUE, show_S = TRUE, hr = TRUE, tick = TRUE,
    theme, mar = c(5, 4, 4, 4) + 0.1, bty = "o",
    lower = 0, upper = 1, pos = 0, horizontal=TRUE, ...)

\method{print}{opticut1}(x, cut, sort, digits, ...)
\method{print}{opticut}(x, digits, ...)
\method{print}{summary.opticut}(x, cut, sort, digits, ...)
\method{summary}{opticut}(object, ...)

\method{as.data.frame}{opticut}(x,
    row.names = NULL, optional = FALSE, cut, sort, ...)
\method{as.data.frame}{summary.opticut}(x,
    row.names = NULL, optional = FALSE, cut, sort, ...)
}
\arguments{
  \item{formula}{
two sided model formula, response species data (matrix,
or possible a vector for single species case) in the left-hand side,
model terms for modifying effects in the right-hand side
(its structure depending on the underlying functions).
For example, in the most basic Gaussian case it can be
\code{y ~ 1} (mo modifying variables) or \code{y ~ x}
(with modifying variables).
Centering the modifying terms (or choosing the origin wisely)
is generally recommended (especially for Gaussian distribution
where linear predictors are additive on the response scale)
because the relative abundance contrast is estimated at the origin (0).
}
  \item{data}{
an optional data frame, list or environment containing the variables
in the model. If not found in data, the variables are taken from
\code{parent.frame()}, typically the environment
from which \code{opticut} is called.
}
  \item{strata}{
vector (usually a factor), unique values define partitions.
It can also be a matrix with rows as observations and
binary partitions as columns.
}
  \item{dist}{
character or function, a distribution to fit.
If character, it can follow one of these patterns: \code{"family"},
or \code{"family:link"} when appropriate (there is a \code{link}
argument in the underlying function, or the link
can be specified via the \code{family} argument).
See Details and Examples.
}
  \item{comb}{
character, how to define the binary partitions.
\code{"rank"} uses \code{\link{rankComb}},
\code{"all"} uses \code{\link{allComb}}.
}
  \item{sset}{
an optional vector specifying a subset of observations to be used in the fitting process. \code{NULL} means no subset taken.
}
  \item{cl}{
a cluster object, or an integer for multiple cores in parallel computations
(this does nothing on Windows).
}
  \item{Y}{
numeric vector of observations for \code{opticut1},
vector or community matrix for \code{opticut.default}.
}
  \item{X}{
numeric, design matrix. Can be missing, in which case an intercept-only model
is assumed.
}
  \item{Z}{
factor (must have at least 2 unique levels,
this triggers \code{\link{rankComb}}),
or a design matrix (custom matrix or as returned by \code{\link{allComb}}.
}
  \item{x, object}{
object to plot, print, summarize. For \code{fix_levels} it needs to be a factor.
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
\code{TRUE}, wile any other numeric value is equivalent to \code{FALSE}).
}
  \item{show_I}{
logical, if indicator values should be shown.
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
  \item{ylab, xlab, las, ylim}{
graphical arguments, see \code{\link{plot}}.
By default, \code{las} is 1 when \code{horizontal = TRUE} and
2 when \code{horizontal = FALSE}.
}
  \item{bty}{
Character, determines the type of box which is drawn about plots,
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
or or \code{NULL} (all species, default).
}
  \item{sep}{
a character string to separate the sub-strings in factor levels.
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
  \item{pos_only}{
logical, best partition normally returns the original variable without
recognizing the direction of the association.
\code{pos_only = TRUE} returns values where negative associations are
taken into account and \code{1} indicates strata of positive association.
This is only important when \code{comb} is not \code{"rank"}.
}
  \item{gnew, xnew}{
new values for \code{strata} and modifiers (right-hand-side of formula)
to predict for, or \code{NULL}. If not \code{NULL}, both must be provided
(available for \code{comb = "rank"} models only).
}
  \item{\dots}{
other arguments passed to the underlying functions.
}
}
\details{
Currently available distributions:

\describe{
\item{\code{"gaussian"}}{real valued continuous observations, e.g. biomass,
uses \code{\link[stats]{lm}} of the stats package.}

\item{\code{"poisson"}}{Poisson count data,
uses \code{\link[stats]{glm}} of the stats package.}

\item{\code{"binomial"}}{presence-absence (detection-nondetection) type data,
uses \code{\link[stats]{glm}} of the stats package.}

\item{\code{"negbin"}}{overdispersed Negative Binomial count data,
uses \code{\link[MASS]{glm.nb}} of the MASS package.}

\item{\code{"beta"}}{continuous response in the unit interval, e.g. percent cover,
uses \code{\link[betareg]{betareg}}  of the betareg package.}

\item{\code{"zip"}}{zero-inflated Poisson counts,
indicative properties are tested as part of the abundance model,
uses \code{\link[pscl]{zeroinfl}} of the pscl package.
The zero-inflation component refers to the probability of 0.}

\item{\code{"zinb"}}{zero-inflated Negative Binomial counts,
indicative properties are tested as part of the abundance model,
uses \code{\link[pscl]{zeroinfl}} of the pscl package.
The zero-inflation component refers to the probability of 0.}

\item{\code{"zip2"}}{zero-inflated Poisson counts,
indicative properties are tested as part of the zero-model,
uses \code{\link[pscl]{zeroinfl}} of the pscl package.
The zero-inflation component refers to the probability of 1
to be consistent with other methods regarding positive and negative effects.}

\item{\code{"zinb2"}}{zero-inflated Negative Binomial counts,
indicative properties are tested as part of the zero-model,
uses \code{\link[pscl]{zeroinfl}} of the pscl package.
The zero-inflation component refers to the probability of 1
to be consistent with other methods regarding positive and negative effects.}

\item{\code{"ordered"}}{response measured on ordinal scale,
e.g. ordinal vegetation cover,
uses \code{\link[MASS]{polr}} of the MASS package.}

\item{\code{"rsf"}}{presence-only data using resource selection
functions (RSF) as explained in  \code{\link[ResourceSelection]{rsf}}
in the ResourceSelection package, assuming global availability (\code{m = 0}).
The \code{"rsf"} works only for single species using \code{opticut1}
because 'presence-only' type data cannot be kept in a single
matrix-like object for multiple species.
Intercept only model is accepted for \code{"rsf"}.}

\item{\code{"rspf"}}{presence-only data using resource selection
probability functions (RSPF)
as explained in \code{\link[ResourceSelection]{rspf}}
in the ResourceSelection package, assuming global availability (\code{m = 0}).
The \code{"rspf"} works only for single species using \code{opticut1}
because 'presence-only' type data cannot be kept in a single
matrix-like object for multiple species.
Intercept only model is not accepted for \code{"rspf"}.}
}

Custom distributions can be defined, see Examples.

\code{fix_levels} is a utility function for replacing characters in
factor levels that are identical to the value of the
\code{getOption("ocoptions")$collapse} value.
This case can lead to an error when specifying the \code{strata} argument,
and the \code{fix_levels} can help.
}
\section{Warning}{
The use of the \code{opticut1} function is generally discouraged:
some of the internal checks are not guaranteed to
flag issues when the formula-to-model-matrix translation is side-stepped
(this is what is happening when the modifier variables are supplied
as \code{X} argument in \code{opticut1}).
Use the \code{opticut} with a single species instead.
}
\value{
\code{opticut1} returns an object of class opticut1, that is a modified
data frame with additional attributes.

\code{opticut} returns an object of class opticut, that is a list
with the following components:
\describe{
\item{\code{"call"}}{the function call.}
\item{\code{"species"}}{a list of species specific opticut1 objects.}
\item{\code{"X"}}{modifying variables as model matrix.}
\item{\code{"Y"}}{response, single species vector or matrix.}
\item{\code{"strata"}}{defines the partitions.}
\item{\code{"nobs"}}{sample size.}
\item{\code{"sset"}}{subset, if specified.}
\item{\code{"nsplit"}}{number of binary splits considered.}
\item{\code{"dist"}}{distribution.}
\item{\code{"comb"}}{combination type.}
\item{\code{"failed"}}{IDs for failed species models.}
\item{\code{"collapse"}}{character used for combining partition labels.}
}

\code{fix_levels} returns a factor with modified levels.

The \code{strata} method extracts the \code{strata} argument
as factor. The method finds unique row combinations
when custom matrix is supplied for \code{strata}.

The \code{print} and \code{summary} methods are called for their side effects.
The summary shows the following information:
best supported split, strength and sign of association,
indicator value (I), expected values (mu0, mu1), log likelihood ratio (logLR),
and model weights(w).

The \code{plot} method presents the contrasts by species and strata.

The \code{wplot} (weight plot) shows model weights for partitions.

\code{\link{bestpart}} returns a matrix with the best supported
partitions for each species (species as columns).

\code{\link{bestmodel}} returns the best supported model for further
manipulation (e.g. prediction). Note: custom distribution
functions are designed to return only point estimates,
thus the best model cannot be returned. In this case,
use the best partition returned by \code{bestpart} to refit the model.
\code{getMLE} returns a named list corresponding to the best supported
model. The list has the following elements:
\code{coef} is the Maximum Likelihood Estimate (MLE),
\code{vcov} is the variance-covariance matrix for the MLE,
\code{dist} is the distribution inherited from input \code{object}.

\code{\link{fitted}} returns expected values on the predictor scale
for the observations as a matrix (number of observations by number of species).
\code{\link{predict}} returns \code{fitted} values when both \code{gnew}
and \code{xnew} are \code{NULL}, or corresponding point predictions
(expected values) on the predictor scale
(available for \code{comb = "rank"} models only).

The coercion method \code{\link{as.data.frame}} return a data frame.
}
\references{
Kemencei, Z., Farkas, R., Pall-Gergely, B., Vilisics, F., Nagy, A., Hornung,
E. & Solymos, P. (2014): Microhabitat associations of land snails in
forested dolinas: implications for coarse filter conservation.
Community Ecology 15:180--186.
<doi:10.1556/ComEc.15.2014.2.6>
}
\author{
Peter Solymos <solymos@ualberta.ca>
}
\seealso{
\code{\link{allComb}}, and \code{\link{rankComb}}
for partitioning algorithms.

\code{\link{beta2i}} for indicator value (I) calculations in summaries.

\code{\link{bestmodel}}, \code{\link{bestpart}}, and \code{\link{uncertainty}}
for manipulating fitted objects.

\code{\link{ocoptions}} on how to set some of the global options
related to the presentation of the results in the package.

\code{\link{optilevels}} for finding the optimal number of factor levels.
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
summary(opticut(spp, strata=gr, dist="gaussian", comb="all")) # alternative

## rank based partitions
summary(ocrank <- opticut(spp ~ 1, strata=gr, dist="gaussian", comb="rank"))
summary(opticut(spp, strata=gr, dist="gaussian", comb="rank")) # alternative

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
(m1 <- opticut(Y ~ x2, strata=x0, dist="poisson", comb="all"))
summary(m1)
## show all species
summary(m1, cut=0)
## plot best partitions and indicator values
plot(m1)
## model weights for all species
wplot(m1)
## different ways of plotting weights for single species
wplot(m1$species[[1]])
wplot(m1, which = 1)

## rank based partitions, multiple species
summary(m2 <- opticut(Y ~ x2, strata=x0, dist="poisson", comb="rank"))

## best partition
head(bestpart(m2))

## best model
mods <- bestmodel(m2)
mods
## explore further
confint(mods[[1]])

## MLE and variance-covariance matrix (species 1)
getMLE(m2, which = 1)

## fitted values
head(fitted(m2))
## prediction for new data
head(predict(m2, gnew=x0, xnex=data.frame(x2=x2)))

\dontrun{
## --- Zero-inflated Negative Binomial
## dolina example
data(dolina)
## stratum as ordinal
dolina$samp$stratum <- as.integer(dolina$samp$stratum)
## filter species to speed up things a bit
Y <- dolina$xtab[,colSums(dolina$xtab > 0) >= 20]
## opticut results, note the cloglog link function
dol <- opticut(Y ~ stratum + lmoist + method, data=dolina$samp,
    strata=dolina$samp$mhab, dist="zinb:cloglog")
summary(dol)
## vertical plot orientation
plot(dol, horizontal=FALSE, pos=1, upper=0.8)

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
summary(dol2)
}

## current collapse value
getOption("ocoptions")$collapse
## factor levels sometimes need to be manipulated
## before feeding it to opticut
fix_levels(as.factor(c("A b", "C d")), sep=":")
fix_levels(as.factor(c("A b", "C d")), sep="")
}
\keyword{ models }
\keyword{ utilities }
