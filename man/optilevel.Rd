\name{optilevel}
\alias{optilevel}
\title{
Optimal number of factor levels.
}
\description{
Finds the optimal number of factor levels given the data and a model.
}
\usage{
optilevel(y, x, alpha = 0, dist = "gaussian", ...)
}
\arguments{
  \item{y}{
vector of observations.
}
  \item{x}{
a factor or a matrix of proportions (i.e. the values 0 and 1 should have
consistent meaning across the columns, often through a unit sum constraint).
}
  \item{alpha}{
weighting factor for calculating information criteria for model selection
(i.e. IC = (1-alpha)*AIC + alpha*BIC, also referred to as CAIC: consistent AIC).
}
  \item{dist}{
character, distribution argument passed to underlying functions,
see listed on the help page of \code{\link{opticut}}
(except for \code{dist = "rsf"} and \code{dist = "rspf"}
which does not have a sound theoretical
basis due to non-identifiability of piece-wise constant RSPF models).
}
  \item{\dots}{
other arguments passed to the underlying functions, see \code{\link{opticut}}.
}
}
\value{
A list with the following elements:

\code{delta}: delta IC values along the selection path considering best models.

\code{ic}: IC values along the selection path considering best models.

\code{coef}: matrix of coefficients (linear predictor scale)
along the selection path considering best models.

\code{rank}: matrix ranks based on the coefficients
along the selection path considering best models.
Ranking uses the default \code{ties.method = "average"} in \code{\link{rank}}.

\code{deltalist}: delta IC values along the selection path
considering all competing models.

\code{iclist}: IC values along the selection path
considering all competing models.

\code{coeflist}: matrix of coefficients (linear predictor scale)
along the selection path considering all competing models.

\code{ranklist}: matrix ranks based on the coefficients
along the selection path considering all competing models.

\code{levels}: list of (merged) factor levels along the selection path
considering best models.
}
\author{
Peter Solymos <solymos@ualberta.ca>
}
\seealso{
\code{\link{opticut}} for finding best binary partitions.

\code{\link{wrsi}} for weighted relative suitability index.
}
\examples{
## --- Factor levels with Gaussian distribution
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

ol <- optilevel(spp[,"Species3"], gr)
ol[c("delta", "coef", "rank", "levels")]

## get the final factor level
gr1 <- gr
levels(gr1) <- ol$level[[length(ol$level)]]
table(gr, gr1)

## compare the models
o0 <- lm(spp[,"Species3"] ~ gr - 1)
o1 <- lm(spp[,"Species3"] ~ gr1 - 1)
data.frame(AIC(o0, o1), delta=AIC(o0, o1)$AIC - AIC(o0))
ol$delta # should be identical

## --- Proportions with Poisson distribution
## simulation
set.seed(123)
n <- 500 # number of observations
k <- 5 # number of habitat types
b <- c(-1, -0.2, -0.2, 0.5, 1)
names(b) <- LETTERS[1:k]
x <- replicate(k, exp(rnorm(n)))
x <- x / rowSums(x) # proportions
X <- model.matrix(~.-1, data=data.frame(x))
lam <- exp(drop(X \%*\% b))
y <- rpois(n, lam)

z <- optilevel(y, x, dist="poisson")

## estimates
plogis(z$coef)
plogis(b)
## optimal classification
z$rank

## get the final matrix
x1 <- mefa4::groupSums(x, 2, z$levels[[length(z$levels)]])
head(x)
head(x1)

## compare the models
m0 <- glm(y ~ x - 1, family="poisson")
m1 <- glm(y ~ x1 - 1, family="poisson")
data.frame(AIC(m0, m1), delta=AIC(m0, m1)$AIC - AIC(m0))
z$delta # should be identical
}
\keyword{ manip }
\keyword{ models }
