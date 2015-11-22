\name{opticut}
\alias{opticut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
opticut(formula, data, strata, dist = "gaussian", comb = c("rank", "all"), cl = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{strata}{
%%     ~~Describe \code{strata} here~~
}
  \item{dist}{
%%     ~~Describe \code{dist} here~~
}
  \item{comb}{
%%     ~~Describe \code{comb} here~~
}
  \item{cl}{
%%     ~~Describe \code{cl} here~~
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
function (formula, data, strata, dist = "gaussian", comb = c("rank", 
    "all"), cl = NULL, ...) 
{
    comb <- match.arg(comb)
    if (missing(data)) 
        data <- parent.frame()
    mf <- match.call(expand.dots = FALSE)
    mm <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, mm)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    Y <- model.response(mf, "numeric")
    Y <- data.matrix(Y)
    if (is.null(colnames(Y))) 
        colnames(Y) <- paste("Species", seq_len(ncol(Y)))
    if (any(duplicated(colnames(Y)))) 
        stop("Duplicate colnames found in LHS")
    ff <- formula
    ff[[2]] <- NULL
    mt <- terms(ff, data = data)
    X <- model.matrix(mt, mf)
    if (missing(strata)) 
        stop("strata is missing")
    if (is.null(dim(strata))) {
        if (comb == "rank") {
            Z <- droplevels(as.factor(strata))
            levels(Z) <- gsub(" +", "", levels(Z))
        }
        if (comb == "all") {
            Z <- allComb(strata)
        }
    }
    else {
        Z <- as.matrix(strata)
        comb <- NA
    }
    Y <- data.matrix(Y)
    if (!is.function(dist)) 
        if (dist == "rspf" && ncol(Y) > 1L) 
            stop("rspf is only available for single species in RHS")
    if (is.null(cl)) {
        if (ncol(Y) > 1L && interactive() && require(pbapply)) {
            res <- pbapply(Y, 2, function(yy, ...) opticut1(Y = yy, 
                X = X, Z = Z, dist = dist, ...), ...)
        }
        else {
            res <- apply(Y, 2, function(yy, ...) opticut1(Y = yy, 
                X = X, Z = Z, dist = dist, ...), ...)
        }
    }
    else {
        library(parallel)
        if (inherits(cl, "cluster")) {
            clusterExport(cl, c("opticut1", ".opticut1", "checkComb", 
                "allComb", "kComb", "rankComb", "oComb"))
            e <- new.env()
            assign("dist", dist, envir = e)
            assign("X", X, envir = e)
            assign("Z", X, envir = e)
            clusterExport(cl, c("X", "Z", "dist"), envir = e)
            res <- parApply(cl, Y, 2, function(yy, ...) opticut1(Y = yy, 
                X = X, Z = Z, dist = dist, ...), ...)
            clusterEvalQ(cl, rm(list = c("opticut1", ".opticut1", 
                "X", "Z", "dist", "checkComb", "allComb", "kComb", 
                "rankComb", "oComb")))
        }
        else {
            if (cl < 2) 
                stop("cl must be at least 2 for forking")
            res <- mclapply(1:ncol(Y), function(i, ...) opticut1(Y = Y[, 
                i], X = X, Z = Z, dist = dist, ...), ...)
        }
    }
    out <- list(call = match.call(), species = res, X = X, Y = Y, 
        strata = Z, dist = dist, comb = comb)
    class(out) <- "opticut"
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
