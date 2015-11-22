\name{summary.opticut}
\alias{summary.opticut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.opticut(object, cut = 2, sort = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{cut}{
%%     ~~Describe \code{cut} here~~
}
  \item{sort}{
%%     ~~Describe \code{sort} here~~
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
function (object, cut = 2, sort = TRUE, ...) 
{
    spp <- lapply(object$species, function(z) as.matrix(z[order(z$w, 
        decreasing = TRUE)[1L], ]))
    sppmat <- t(sapply(spp, function(z) as.matrix(z)))
    hab <- sapply(spp, rownames)
    colnames(sppmat) <- colnames(object$species[[1L]])
    res <- data.frame(split = hab, sppmat)
    res$assoc <- parseAssoc(res)
    if (sort) 
        res <- res[order(res$split, 1 - res$I, decreasing = FALSE), 
            ]
    res <- res[res$logLR >= min(max(res$logLR), cut), , drop = FALSE]
    res$logL <- NULL
    object$summary <- res
    object$nsplit <- if (is.factor(object$strata)) 
        nlevels(object$strata)
    else ncol(object$strata)
    object$missing <- length(object$species) - nrow(res)
    class(object) <- c("summary.opticut")
    object
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
