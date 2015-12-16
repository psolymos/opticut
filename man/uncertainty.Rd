\name{uncertainty}
\alias{uncertainty}
\alias{uncertainty.opticut}
\alias{print.uncertainty1}
\alias{plot.uncertainty1}
\alias{wplot.uncertainty1}
\alias{print.uncertainty}
\alias{print.summary.uncertainty}
\alias{summary.uncertainty}
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

\method{print}{uncertainty1}(x, ...)
\method{plot}{uncertainty1}(x, ...)
\method{wplot}{uncertainty1}(x, ylab, ...)
\method{print}{uncertainty}(x, ...)
\method{print}{summary.uncertainty}(x, ...)
\method{summary}{uncertainty}(object, level = 0.95, ...)
}
\arguments{
  \item{object}{
fitted model object from \code{\link{opticut}} (or an output from
\code{uncertainty} for the \code{summary} method).
}
  \item{which}{
numeric, character (can be a vector) defining
a subset of species from the fitted object,
or or \code{NULL} (all species, default).
}
  \item{type}{
character, describing the type of uncertainty calculation for
indicator value \code{I}, and expected values \code{mu0}, and \code{mu1}
(\code{I = 1 - min(mu0, mu1) / max(mu0, mu1)}).
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
an object to be printed or plotted.
}
  \item{ylab}{
label for the y axis.
}
  \item{level}{
the confidence level required.
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
lam <- exp(0.5 + 1*x1 + -0.2*x2)
Y <- rpois(n, lam)
oc <- opticut(Y ~ x2, strata=x0, dist="poisson", comb="rank")

uc1 <- uncertainty(oc, type="asymp", B=1000)
uc2 <- uncertainty(oc, type="boot", B=20)
## use user-supplied indices
B <- replicate(20, sample.int(n, replace=TRUE))
uc3 <- uncertainty(oc, type="multi", B=B)

## selection probabilities (p) for best supported split
## and confidence intervals for I
summary(uc3)

## individual species results
uc3$uncertainty

## density plots
## note: multimodality is due to uncertainty in partitions
plot(uc3$uncertainty[[1]])
## selection frequencies for partitions
## (only for type = 'multi')
wplot(uc3$uncertainty[[1]])
}
\keyword{ manip }
\keyword{ models }
