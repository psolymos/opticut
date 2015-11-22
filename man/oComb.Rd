\name{oComb}
\alias{oComb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
oComb(x, collapse = " ")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{collapse}{
%%     ~~Describe \code{collapse} here~~
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
function (x, collapse = " ") 
{
    if (length(x) < 2L) 
        stop("length of x must be >1")
    if (is.null(names(x))) 
        names(x) <- seq_len(length(x))
    o <- x[order(x, decreasing = FALSE)]
    out <- diag(1L, length(o))
    out[upper.tri(out)] <- 1L
    out <- out[, -ncol(out)]
    rownames(out) <- names(o)
    colnames(out) <- seq_len(ncol(out))
    for (i in seq_len(ncol(out))) {
        colnames(out)[i] <- paste(names(x)[names(x) \%in\% rownames(out)[out[, 
            i] > 0]], collapse = collapse)
    }
    attr(out, "rank") <- o
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
