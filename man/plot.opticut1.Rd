\name{plot.opticut1}
\alias{plot.opticut1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.opticut1(x, cut = 2, ylim = c(-1, 1), ylab = "Model weight * Association", xlab = "Partitions", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{cut}{
%%     ~~Describe \code{cut} here~~
}
  \item{ylim}{
%%     ~~Describe \code{ylim} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
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
function (x, cut = 2, ylim = c(-1, 1), ylab = "Model weight * Association", 
    xlab = "Partitions", ...) 
{
    w <- x$w * x$assoc
    names(w) <- rownames(x)
    if (!any(x$logLR >= cut)) {
        warning("All logLR < cut: cut ignored")
    }
    else {
        w <- w[x$logLR >= cut]
    }
    COL <- c(colorRampPalette(c("red", "yellow"))(10), colorRampPalette(c("yellow", 
        "green"))(10))
    br <- seq(-1, 1, 0.1)
    barplot(rep(0, length(w)), width = 1, space = 0, col = COL[as.integer(cut(w, 
        breaks = seq(-1, 1, 0.1)))], ylim = ylim, xlab = xlab, 
        ylab = ylab, ...)
    lines(rep(which.max(abs(w)) - 0.5, 2), c(-1, 1), col = "grey", 
        lwd = 2)
    barplot(w, width = 1, space = 0, col = COL[as.integer(cut(w, 
        breaks = seq(-1, 1, 0.1)))], ylim = ylim, xlab = "", 
        ylab = "", add = TRUE, ...)
    abline(0, 0)
    box()
    invisible(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
