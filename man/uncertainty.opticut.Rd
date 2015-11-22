\name{uncertainty.opticut}
\alias{uncertainty.opticut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
uncertainty.opticut(object, which = NULL, type = c("asymp", "boot", "multi"), B = 99, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{which}{
%%     ~~Describe \code{which} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
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
function (object, which = NULL, type = c("asymp", "boot", "multi"), 
    B = 99, ...) 
{
    require(pbapply)
    type <- match.arg(type)
    B <- as.integer(B)
    if (B < 1) 
        stop("B must be > 0")
    linkinv <- .opticut1(Y = object$Y[, 1L], X = object$X, Z1 = NULL, 
        dist = object$dist, ...)$linkinv
    m <- .extractOpticut(object, which, boot = FALSE, internal = TRUE, 
        full_model = TRUE, best = TRUE, ...)
    spp <- names(m)
    out <- list()
    if (type == "asymp") {
        for (i in spp) {
            require(MASS)
            k <- which.max(object$species[[i]]$logLR)
            bm <- rownames(object$species[[i]])[k]
            m1 <- m[[i]]
            cf <- mvrnorm(B, coef(m1), vcov(m1))[, c(1L, 2L)]
            cf <- rbind(coef(m1)[c(1L, 2L)], cf)
            cf0 <- linkinv(cf[, 1L])
            cf1 <- linkinv(cf[, 1L] + cf[, 2L])
            I <- 1 - (pmin(cf0, cf1)/pmax(cf0, cf1))
            out[[i]] <- data.frame(best = bm, I = I, mu0 = cf0, 
                mu1 = cf1)
        }
    }
    if (type == "boot") {
        for (i in spp) {
            k <- which.max(object$species[[i]]$logLR)
            bm <- rownames(object$species[[i]])[k]
            cf <- t(pbsapply(seq_len(B), function(z) {
                .extractOpticut(object, i, boot = TRUE, internal = TRUE, 
                  full_model = FALSE, best = TRUE, ...)[[1L]]$coef[c(1L, 
                  2L)]
            }))
            cf <- rbind(coef(m[[i]])[c(1L, 2L)], cf)
            cf0 <- linkinv(cf[, 1L])
            cf1 <- linkinv(cf[, 1L] + cf[, 2L])
            I <- 1 - (pmin(cf0, cf1)/pmax(cf0, cf1))
            out[[i]] <- data.frame(best = bm, I = I, mu0 = cf0, 
                mu1 = cf1)
        }
    }
    if (type == "multi") {
        for (i in spp) {
            if (object$comb == "all") 
                stop("comb='all' is no good, use 'rank' instead")
            k <- which.max(object$species[[i]]$logLR)
            bm <- character(B + 1L)
            bm[1L] <- rownames(object$species[[i]])[k]
            mat <- matrix(NA, B + 1L, 3)
            colnames(mat) <- c("I", "mu0", "mu1")
            tmp <- as.numeric(object$species[[i]][k, -1L])
            names(tmp) <- colnames(object$species[[i]])[-1L]
            mat[1L, ] <- tmp[c("I", "mu0", "mu1")]
            pb <- startpb(0, B)
            on.exit(closepb(pb))
            for (j in seq_len(B)) {
                mod <- .extractOpticut(object, i, boot = TRUE, 
                  internal = FALSE, best = FALSE)[[1L]]
                k <- which.max(mod$logLR)
                bm[j + 1L] <- rownames(mod)[k]
                tmp <- as.numeric(mod[k, -1L])
                names(tmp) <- colnames(mod)[-1L]
                mat[j + 1L, ] <- tmp[c("I", "mu0", "mu1")]
                setpb(pb, j)
            }
            out[[i]] <- data.frame(best = bm, mat)
        }
    }
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
