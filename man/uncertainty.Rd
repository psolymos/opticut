\name{uncertainty}
\alias{uncertainty}
\alias{uncertainty.opticut}
\alias{uncertainty.multicut}
\alias{print.uncertainty1}
\alias{print.uncertainty}
\alias{print.summary.uncertainty}
\alias{summary.uncertainty}
\alias{strata.uncertainty}
\alias{bestpart.uncertainty}
\alias{bestpart.uncertainty1}
\alias{check_strata}
\alias{as.data.frame.uncertainty}
\alias{as.data.frame.summary.uncertainty}
\alias{bsmooth}
\alias{bsmooth.uncertainty1}
\alias{bsmooth.uncertainty}
\alias{subset.uncertainty}
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
\method{uncertainty}{multicut}(object,
    which = NULL, type = c("asymp", "boot"),
    B = 99, cl = NULL, ...)

check_strata(x, mat)
\method{strata}{uncertainty}(object, ...)
\method{subset}{uncertainty}(x, subset=NULL, ...)

\method{bestpart}{uncertainty}(object, ...)
\method{bestpart}{uncertainty1}(object, ...)

\method{print}{uncertainty1}(x, ...)
\method{print}{uncertainty}(x, ...)
\method{print}{summary.uncertainty}(x, sort, digits, ...)
\method{summary}{uncertainty}(object, level = 0.95, ...)

\method{as.data.frame}{uncertainty}(x,
    row.names = NULL, optional = FALSE, sort, ...)
\method{as.data.frame}{summary.uncertainty}(x,
    row.names = NULL, optional = FALSE, sort, ...)

\method{bsmooth}{uncertainty1}(object, ...)
\method{bsmooth}{uncertainty}(object, ...)
}
\arguments{
  \item{object}{
fitted model object from \code{\link{opticut}}
which does not contain extra arguments as part of \code{...}
(or an output from \code{uncertainty} for the \code{summary} method).
}
  \item{which}{
numeric or character (can be a vector) defining
a subset of species from the fitted object,
or or \code{NULL} (all species, default).
}
  \item{type}{
character, describing the type of uncertainty calculation.
See Details.
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
  \item{sort}{
logical value indicating if species
should be meaningfully sorted, the default is \code{TRUE}.
}
  \item{digits}{
numeric, number of significant digits in output.
}
  \item{mat}{
a matrix with resampling indices (rows as samples, columns as iterations).
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
missing values are not accepted.
}
\item{\dots}{
other arguments passed to the underlying functions.
}
}
\details{
Uncertainty is calculated for
indicator value \code{I}, and expected values
(\code{mu0}, and \code{mu1} for opticut, and \code{mu_*} for multicut objects).

\code{"asymp"}: asymptotic distribution is based on best partition
found for the input object (this option
is unavailable for custom distribution functions because
it requires the Hessian matrix).
This type is available for both opticut and multicut objects.

\code{"boot"}: non-parametric bootstrap distribution
based on best partition found for the input object.
This type is available for both opticut and multicut objects.

\code{"multi"}: non-parametric bootstrap distribution
based on best partition found for the bootstrap data (i.e.
the model ranking is re-evaluated each time).
\code{"multi"} works only if \code{comb = "rank"} in the
\code{\link{opticut}} call.
This type is not available for multicut objects.
}
\value{
\code{uncertainty} returns an object of class 'uncertainty'.
The \code{uncertainty} element of the object is a list with species specific
output as elements (object class 'uncertainty1').
Each 'uncertainty1' output is a data frame with 4 columns:
\code{best} partition, indicator value \code{I},
and expected values
(\code{mu0}, and \code{mu1} for opticut, and \code{mu_*} for multicut objects).

\code{check_strata} returns a logical vector checking if
all original strata from the input object are represented
by resampling indices. Number of strata are attached as attributes
for further diagnostics.

The summary method prints the name of the best supported split,
selection frequency (R, reliability), indicator values (I, based on
the distribution of values within the best supported split with highest
reliability) and conficence interval for I (based on \code{level}).

The \code{subset} method subsets the species in the uncertainty object.

\code{\link{bestpart}} finds the selection frequencies for
strata as best partitions (number of strata times number of species).

The coercion method \code{\link{as.data.frame}} returns a data frame.

The \code{bsmooth} method returns bootstrap smoothed results
for each strata (not available for multicut based uncertainty objects).
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
\code{\link{opticut}} and \code{\link{multicut}} for the
user interface of the input objects.
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
(uc1 <- uncertainty(oc, type="asymp", B=999))
summary(uc1)
## bootstrap-based confidence intervals
(uc2 <- uncertainty(oc, type="boot", B=19))
summary(uc2)

## use user-supplied indices
## multi-model bootstrap based uncertainties
B <- replicate(25, sample.int(n, replace=TRUE))
check_strata(oc, B) # check representation
(uc3 <- uncertainty(oc, type="multi", B=B))
summary(uc3)

## best partitions:
## selection frequencies for strata and species
bestpart(uc3)
heatmap(bestpart(uc3), scale="none", col=occolors()(25))

## bootstrap smoothed predictions per strata
bsmooth(uc3)
heatmap(bestpart(uc3), scale="none", col=occolors()(25))

## individual species results
uc3$uncertainty
bestpart(uc3$uncertainty[[1]])
bsmooth(uc3$uncertainty[[1]])

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

## multicut based uncertainty
mc <- multicut(Y ~ x2, strata=x0, dist="poisson")

## asymptotic confidence intervals
(muc1 <- uncertainty(mc, type="asymp", B=999))
summary(muc1)
bestpart(muc1)

## bootstrap-based confidence intervals
(muc2 <- uncertainty(mc, type="boot", B=19))
summary(muc2)
bestpart(muc2)

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

## See how indicator value changes with different partitions
## (and why it is the wrong metric to use in this calse)
with(ucdol$uncertainty[["pvic"]],
    boxplot(I ~ best, col="gold", ylab="Indicator value"))
## What we should calculate is the bootstrap smoothed mean of the
## expected value and its confidence intervals
bs <- bsmooth(ucdol$uncertainty[["pvic"]])
boxplot(t(bs), ylab="Expected value")
cbind(Mean=rowMeans(bs), t(apply(bs, 1, quantile, probs=c(0.025, 0.975))))

## A more interesting simulated example for bootstrap smoothing
## and comparing opticut vs. multicut
set.seed(1)
n <- 2000
x <- sort(runif(n, -8, 8))
p <- plogis(0.5 + -0.1 * x + -0.2 * x^2)
y <- rbinom(n, 1, p)
d <- diff(range(x))/10
br <- seq(min(x), max(x), by=d)
g <- cut(x, br, include.lowest=TRUE)
levels(g) <- LETTERS[1:nlevels(g)]
o <- opticut(y ~ 1, strata=g, dist="binomial")
m <- multicut(y ~ 1, strata=g, dist="binomial")
library(parallel)
cl <- makeCluster(2)
uo <- uncertainty(o, type="multi", B=99, cl=cl)
um <- uncertainty(m, type="boot", B=99, cl=cl)
stopCluster(cl)
## bootstrap average for opticut
bs <- bsmooth(uo$uncertainty[[1]])
stat <- cbind(Mean=rowMeans(bs),
    t(apply(bs, 1, quantile, probs=c(0.025, 0.975))))
## bootstrap average for multicut
bsm <- as.matrix(um$uncertainty[[1]][,-(1:2)])
statm <- cbind(Mean=colMeans(bsm),
    t(apply(bsm, 2, quantile, probs=c(0.025, 0.975))))

op <- par(mfrow=c(2,1))
plot(p ~ x, type="l", ylim=c(0,1), main="Binary partitions (opticut)")
abline(v=br, col="grey", lty=3)
lines(br[-1]-0.5*d, stat[,1], col=4)
lines(br[-1]-0.5*d, stat[,2], col=4, lty=2)
lines(br[-1]-0.5*d, stat[,3], col=4, lty=2)
lines(br[-1]-0.5*d, bs[,1], col=2)
legend("topright", bty="n", lty=c(1,1,2,1), col=c(1,4,4,2),
    legend=c("True response","bsmooth","0.95 CI","Best partition"))

plot(p ~ x, type="l", ylim=c(0,1), main="Multi-level model (multicut)")
abline(v=br, col="grey", lty=3)
lines(br[-1]-0.5*d, statm[,1], col=4)
lines(br[-1]-0.5*d, statm[,2], col=4, lty=2)
lines(br[-1]-0.5*d, statm[,3], col=4, lty=2)
legend("topright", bty="n", lty=c(1,1,2), col=c(1,4,4),
    legend=c("True response","bsmooth","0.95 CI"))
par(op)
}
}
\keyword{ manip }
\keyword{ models }
