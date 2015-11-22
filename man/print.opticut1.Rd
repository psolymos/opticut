\name{print.opticut1}
\alias{print.opticut1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.opticut1(x, cut = 2, sort = TRUE, digits, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{cut}{
%%     ~~Describe \code{cut} here~~
}
  \item{sort}{
%%     ~~Describe \code{sort} here~~
}
  \item{digits}{
%%     ~~Describe \code{digits} here~~
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
function (x, cut = 2, sort = TRUE, digits, ...) 
{
    if (missing(digits)) 
        digits <- max(3L, getOption("digits") - 3L)
    xx <- x
    if (sort) 
        xx <- xx[order(xx$I, decreasing = TRUE), ]
    xx$assoc <- parseAssoc(xx)
    if (any(xx$logLR >= cut)) {
        SHOW <- which(xx$logLR >= cut)
        tmp <- if (length(SHOW) > 1L) 
            "Best supported models"
        else "Best supported model"
        TXT <- paste0(tmp, " with logLR >= ", format(cut, digits = digits), 
            ":")
    }
    else {
        SHOW <- 1L
        TXT <- paste0("Best supported model:")
    }
    xx <- xx[SHOW, , drop = FALSE]
    cat("Univariate opticut results, comb = ", attr(x, "comb"), 
        ", dist = ", attr(x, "dist"), "\nI = ", format(xx[1L, 
            "I"], digits = digits), "; w = ", format(xx[1L, "w"], 
            digits = digits), "; H = ", format(attr(x, "H"), 
            digits = digits), "; logL_null = ", format(attr(x, 
            "logL_null"), digits = digits), "\n\n", TXT, "\n", 
        sep = "")
    print.data.frame(xx, digits = digits, ...)
    DROP <- nrow(x) - nrow(xx)
    if (DROP > 0) {
        cat(nrow(x), " binary ", ifelse(nrow(x) > 1, "splits", 
            "split"), " (", DROP, ifelse(DROP > 1, " models", 
            " model"), " not shown)\n", sep = "")
    }
    else {
        cat(nrow(x), "binary", ifelse(nrow(x) > 1, "splits", 
            "split"), "\n")
    }
    cat("\n")
    invisible(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
