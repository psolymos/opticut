## Y is abundance vector
## X is model matrix for nuisance variables
## Z is factor
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
    n <- table(Z)
    n <- structure(as.numeric(n), names=names(n))
    Z0 <- model.matrix(~Z)
    res0 <- .opticut1(Y, X, Z1=NULL,
        linkinv=TRUE, dist=dist, ...)
    res <- .opticut1(Y, X, Z1=Z0[,-1L,drop=FALSE],
        linkinv=TRUE, dist=dist, ...)
    cf <- res$coef
    gamma <- c(cf[1], cf[1] + cf[2:ncol(Z0)])
    mu <- res$linkinv(gamma)
    names(mu) <- levels(Z)
    K <- nlevels(Z)
    fix_fitted <- getOption("ocoptions")$fix_fitted
    mufix <- if (fix_fitted && any(mu < 0))
        mu + abs(min(mu)) else mu
    l <- lorenz(mufix, n)
    bp <- structure(numeric(K), names=names(mu))
    bp[rownames(l)[-seq_len(which.max(l[,"p"] - l[,"L"]))]] <- 1
    bpvec <- bp[as.integer(Z)]
#    resbp <- .opticut1(Y, X, Z1=bpvec,
#        linkinv=TRUE, dist=dist, ...)
#    cfbp <- resbp$coef
    ll <- res$logLik
    scale <- getOption("ocoptions")$scale
    out <- list(
        null=res0$linkinv(res0$coef[1L]),
        mu=mu,
        bestpart=bp,
#        beta0=cfbp[,1L],
#        beta1=cfbp[,2L],
#        gamma0=cfbp[,1L],
#        gamma1=cfbp[,1L] + cfbp[,2L],
#        I=beta2i(cfbp[,2L], scale=scale),
        #I=max(mulink)-min(mulink),
        I=beta2i(max(gamma) - min(gamma), scale=scale),
        coefficients=cf,
        n=n,
        logL=res$logLik,
        logLR=res$logLik-res0$logLik)
    attr(out, "scale") <- scale
    attr(out, "logL_null") <- res0$logLik
#    attr(out, "logL_bp") <- resbp$logLik
    attr(out, "dist") <- if (is.function(dist))
        deparse(substitute(dist)) else .opticut_dist(dist, make_dist=TRUE)
    class(out) <- c("multicut1")
    out
}
