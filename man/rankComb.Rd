\name{rankComb}
\alias{rankComb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rankComb(Y, X, Z, dist = "gaussian", collapse = " ", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{dist}{
%%     ~~Describe \code{dist} here~~
}
  \item{collapse}{
%%     ~~Describe \code{collapse} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Y, X, Z, dist = "gaussian", collapse = " ", ...) 
{
    if (!is.factor(Z)) 
        stop("Z must be a factor")
    Z0 <- model.matrix(~Z)
    m <- .opticut1(Y, X, Z1 = Z0[, -1L, drop = FALSE], linkinv = TRUE, 
        dist = dist, ...)
    lc <- c(m$coef[1], m$coef[1] + m$coef[2:ncol(Z0)])
    names(lc) <- levels(Z)
    x <- rank(-lc)
    oc <- oComb(x, collapse = collapse)
    out <- oc[match(Z, rownames(oc)), ]
    attr(out, "est") <- m$linkinv(lc)
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
