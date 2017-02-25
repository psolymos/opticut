## Y is abundance vector
## X is model matrix for nuisance variables
## Z is design matrix for binary splits or a factor (using rankComb)
multicut1 <-
function(Y, X, Z, dist="gaussian", sset=NULL, ...)
{
    if (missing(X))
        X <- matrix(1L, length(Y), 1L)
    X <- data.matrix(X)
    if (is.null(rownames(X)))
        rownames(X) <- seq_len(nrow(X))
    if (!is.factor(Z))
        stop("Z must be a factor")
    if (!is.null(sset)) {
        Y <- Y[sset]
        X <- X[sset,,drop=FALSE]
        Z <- Z[sset]
    }
    Z0 <- model.matrix(~Z)
    res0 <- .opticut1(Y, X, Z1=NULL,
        linkinv=TRUE, dist=dist, ...)
    res <- .opticut1(Y, X, Z1=Z0[,-1L,drop=FALSE],
        linkinv=TRUE, dist=dist, ...)
    cf <- res$coef
    mu <- res$linkinv(c(cf[1], cf[1] + cf[2:ncol(Z0)]))
    names(mu) <- levels(Z)
    ll <- res$logLik
    cfnull <- res0$linkinv(res0$coef[1L])
    out <- list(
        null=res0$linkinv(res0$coef[1L]),
        mu=mu,
        coefficients=cf,
        logL=ll,
        logLR=ll-res0$logLik)
    attr(out, "logL_null") <- res0$logLik
    attr(out, "dist") <- if (is.function(dist))
        deparse(substitute(dist)) else dist
    class(out) <- c("multicut1")
    out
}


#set.seed(1234)
#n <- 200
#x0 <- sample(1:4, n, TRUE)
#x1 <- ifelse(x0 %in% 1:2, 1, 0)
#x2 <- rnorm(n, 0.5, 1)
#lam <- exp(0.5 + 0.5*x1 + -0.2*x2)
#Y1 <- rpois(n, lam)
#Y2 <- rpois(n, rev(lam))
#Y <- cbind(Y1, Y2)
#.opticut1=opticut:::.opticut1
#rc <- multicut1(Y1, model.matrix(~x2), as.factor(x0), dist="poisson")
#rc <- multicut(Y, model.matrix(~x2), as.factor(x0), dist="poisson")
#rc <- multicut(Y~x2, strata=as.factor(x0), dist="poisson")
