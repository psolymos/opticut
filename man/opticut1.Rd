\name{opticut1}
\alias{opticut1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
opticut1(Y, X, Z, dist = "gaussian", ...)
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
function (Y, X, Z, dist = "gaussian", ...) 
{
    X <- data.matrix(X)
    if (is.null(rownames(X))) 
        rownames(X) <- seq_len(nrow(X))
    if (is.factor(Z)) {
        Z <- rankComb(Y, X, Z, dist = dist, ...)
        Est <- attr(Z, "est")
        Comb <- "rank"
    }
    else {
        Est <- NA
        Comb <- "all"
    }
    Z <- data.matrix(Z)
    if (is.null(colnames(Z))) 
        colnames(Z) <- paste0("split.", seq_len(ncol(Z)))
    if (!checkComb(Z)) 
        stop("complementary design variables found")
    if (length(unique(c(length(Y), nrow(X), nrow(Z)))) > 1) 
        stop("dimension mismatch")
    N <- ncol(Z)
    res0 <- .opticut1(Y, X, Z1 = NULL, dist = dist, ...)
    cf <- matrix(0, N, length(res0$coef) + 1)
    rownames(cf) <- colnames(Z)
    ll <- numeric(N)
    names(ll) <- colnames(Z)
    for (i in seq_len(N)) {
        res <- .opticut1(Y, X, Z1 = Z[, i], dist = dist, ...)
        cf[i, ] <- res$coef
        ll[i] <- res$logLik
    }
    dll <- ll - max(ll)
    w <- exp(dll)/sum(exp(dll))
    cfnull <- res0$linkinv(res0$coef[1L])
    cf0 <- res0$linkinv(cf[, 1L])
    cf1 <- res0$linkinv(cf[, 1L] + cf[, 2L])
    h <- sign(cf[, 2L])
    I <- 1 - (pmin(cf0, cf1)/pmax(cf0, cf1))
    out <- data.frame(assoc = h, I = I, null = cfnull, mu0 = cf0, 
        mu1 = cf1, logL = ll, logLR = ll - res0$logLik, w = w)
    rownames(out) <- colnames(Z)
    attr(out, "logL_null") <- res0$logLik
    attr(out, "H") <- sum(w^2)
    attr(out, "dist") <- if (is.function(dist)) 
        deparse(substitute(dist))
    else dist
    attr(out, "comb") <- Comb
    attr(out, "est") <- Est
    class(out) <- c("opticut1", "data.frame")
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
