## Y is abundance vector
## X is model matrix for nuisance variables
## Z is factor
##
## does not store .lc_cut1 output because it might fail
## leva that for downstream functions but store ingredients
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
    mulink <- c(cf[1], cf[1] + cf[2:ncol(Z0)])
    mu <- res$linkinv(mulink)
    names(mu) <- levels(Z)
    ll <- res$logLik
    scale <- getOption("ocoptions")$scale
    out <- list(
        null=res0$linkinv(res0$coef[1L]),
        mu=mu,
        #I=max(mulink)-min(mulink),
        I=beta2i(max(mulink) - min(mulink), scale=scale),
        coefficients=cf,
        n=n,
        logL=res$logLik,
        logLR=res$logLik-res0$logLik)
    attr(out, "scale") <- scale
    attr(out, "logL_null") <- res0$logLik
    attr(out, "dist") <- if (is.function(dist))
        deparse(substitute(dist)) else .opticut_dist(dist, make_dist=TRUE)
    class(out) <- c("multicut1")
    out
}
