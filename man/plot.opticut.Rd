\name{plot.opticut}
\alias{plot.opticut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.opticut(x, which = NULL, cut = 2, sort = TRUE, coverage = 0.95, las = 1, ylab = "Model weight * Association", xlab = "Partitions", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{which}{
%%     ~~Describe \code{which} here~~
}
  \item{cut}{
%%     ~~Describe \code{cut} here~~
}
  \item{sort}{
%%     ~~Describe \code{sort} here~~
}
  \item{coverage}{
%%     ~~Describe \code{coverage} here~~
}
  \item{las}{
%%     ~~Describe \code{las} here~~
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
function (x, which = NULL, cut = 2, sort = TRUE, coverage = 0.95, 
    las = 1, ylab = "Model weight * Association", xlab = "Partitions", 
    ...) 
{
    if (!is.null(which) && length(which) == 1L) {
        plot(x$species[[which]], cut = cut, ylab = ylab, xlab = xlab, 
            ...)
    }
    else {
        if (x$comb == "rank") 
            stop("Plot single species (comb='rank'): specify 'which' argument")
        if (!is.null(which) && length(which) > 1L) 
            x$species <- x$species[which]
        COL <- c(colorRampPalette(c("red", "yellow"))(10), colorRampPalette(c("yellow", 
            "green"))(10))
        br <- seq(-1, 1, 0.1)
        xx <- summary(x, cut = cut, sort = sort)
        nsplit <- xx$nsplit
        nspp <- nrow(xx$summary)
        ww <- sapply(xx$species[rownames(xx$summary)], "[[", 
            "w")
        ss <- sapply(xx$species[rownames(xx$summary)], "[[", 
            "assoc")
        ss[ss == 0] <- 1
        ww <- ww * ss
        rownames(ww) <- rownames(xx$species[[1]])
        colnames(ww) <- rownames(xx$summary)
        llr <- sapply(xx$species[rownames(xx$summary)], "[[", 
            "logLR")
        ww[llr < cut] <- 0
        ww <- ww[rowSums(ww) != 0, , drop = FALSE]
        op <- par(las = las)
        plot(0, xlim = c(0, nrow(ww)), ylim = c(ncol(ww), 0), 
            type = "n", axes = FALSE, ann = FALSE, ...)
        title(ylab = ylab, xlab = xlab)
        axis(2, at = 1:ncol(ww) - 0.5, labels = colnames(ww), 
            tick = TRUE)
        axis(1, at = 1:nrow(ww) - 0.5, labels = rownames(ww), 
            tick = TRUE)
        abline(h = 1:ncol(ww) - 0.5)
        abline(v = 0:nrow(ww), col = "lightgrey")
        for (i in 1:ncol(ww)) {
            lines(rep(which.max(abs(ww[, i])), 2) - 0.5, c(-0.5, 
                0.5) + i - 0.5, col = "grey", lwd = 2)
            for (j in 1:nrow(ww)) {
                h <- -ww[j, i] * 0.45
                polygon(c(0, 1, 1, 0) + j - 1, c(0, 0, h, h) + 
                  i - 0.5, col = COL[as.integer(cut(-h, breaks = seq(-1, 
                  1, 0.1)))])
            }
        }
        box()
        par(op)
        invisible(ww)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
