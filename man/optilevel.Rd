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
consistent meaning actoss the columns, often through a unit sum constraint).
}
  \item{alpha}{
weighting factor for calculating information criteria for model selection
(i.e. IC = (1-alpha)*AIC + alpha*BIC, also referred to as CAIC: consistent AIC).
}
  \item{dist}{
character, distribution argument passed to underlying functions,
see listed on the help page of \code{\link{opticut}}.
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
}
\examples{

## show how to get the final merging

## factor and prop based example

## simul (prop) & real (dolina)???

}
\keyword{ manip }
\keyword{ models }
