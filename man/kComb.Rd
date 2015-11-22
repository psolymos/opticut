\name{kComb}
\alias{kComb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
kComb(k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{
%%     ~~Describe \code{k} here~~
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
function (k) 
{
    k <- as.integer(k)
    if (k < 2) 
        stop("k must be at least 2")
    kmax <- floor(k/2)
    s <- seq_len(k)
    clist <- lapply(seq_len(kmax), function(kk) combn(k, kk))
    if (kmax == k/2) {
        COL <- seq_len(ncol(clist[[kmax]])/2)
        clist[[kmax]] <- clist[[kmax]][, COL, drop = FALSE]
    }
    m <- sapply(clist, ncol)
    out <- matrix(0L, k, sum(m))
    z <- 1
    for (i in seq_len(length(clist))) {
        for (j in seq_len(m[i])) {
            out[s \%in\% clist[[i]][, j], z] <- 1L
            z <- z + 1
        }
    }
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
