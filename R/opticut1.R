## Y is abundance vector
## X is model matrix for nuisance variables
## Z is design matrix for binary splits or a factor (using rankComb)
opticut1 <-
function(Y, X, Z, dist="gaussian", sset=NULL, ...)
{
    if (missing(X))
        X <- matrix(1L, length(Y), 1L)
    X <- data.matrix(X)
    if (is.null(rownames(X)))
        rownames(X) <- seq_len(nrow(X))
    if (is.factor(Z)) {
        Z <- rankComb(Y, X, Z, dist=dist, ...)
        Est <- attr(Z, "est")
        Comb <- "rank"
    } else {
            stop("Pssst ... Z must have 0 and 1 values only!")
        Est <- NA
        Comb <- attr(Z, "comb")
        if (is.null(Comb))
            Comb <- NA
    }
    Z <- data.matrix(Z)
    if (is.null(colnames(Z)))
        colnames(Z) <- paste0("split.", seq_len(ncol(Z)))
    if (getOption("ocoptions")$check_comb && !checkComb(Z))
        stop("Guess what! Complementary design variables found:\nuse 'checkComb'")
    if (length(unique(c(length(Y), nrow(X), nrow(Z)))) > 1)
        stop("Khm ... dimension mismatch: check you input.")
    if (is.null(rownames(Z))) {
        #warning("Row names added to binary split matrix Z (it was NULL).\nYou are welcome.")
        rownames(Z) <- apply(Z, 1, paste, collapse="")
    }
    N <- ncol(Z)
    if (!is.null(sset)) {
        Y <- Y[sset]
        X <- X[sset,,drop=FALSE]
        Z <- Z[sset,,drop=FALSE]
    }
    res0 <- .opticut1(Y, X, Z1=NULL, dist=dist, ...)
    cf <- matrix(0, N, length(res0$coef) + 1L)
    rownames(cf) <- colnames(Z)
    ll <- numeric(N)
    names(ll) <- colnames(Z)
    for (i in seq_len(N)) {
        res <- .opticut1(Y, X, Z1=Z[,i], dist=dist, ...)
        cf[i,] <- res$coef
        ll[i] <- res$logLik
    }
    dll <- ll - max(ll)
    w <- exp(dll) / sum(exp(dll))
    cfnull <- res0$linkinv(res0$coef[1L])
    cf0 <- res0$linkinv(cf[,1L])
    cf1 <- res0$linkinv(cf[,1L] + cf[,2L])
    h <- sign(cf[,2L])

    ## AIC weight has a penalty dependent 'midpoint'
    ## only 1 df difference, this length of coef is not important
    #ic <- cbind(ic = -2*ll + getOption("ocoptions")$penalty,
    #    ic_null = -2*res0$logLik)
    #Delta <- t(apply(ic, 1, function(z) z - min(z)))
    #W <- apply(Delta, 1, function(z) exp(-0.5*z[1]) / sum(exp(-0.5*z)))
    ## delta IC with inverse Fisher transform is more intuitive
    ## 0 when there is no support for a partition
    ## 1 when logLR is HUGE
    ## this is trasformed (AIC_null - AIC_m)
    #W <- pmax(0, tanh(2*ll - 2*res0$logLik - getOption("ocoptions")$penalty))

    ## problem with prediction based I is that
    ## - it depends on covariates if not centered
    ## - ordered intercept can be ill-defined in some cases
#    I <- 1 - (pmin(cf0, cf1) / pmax(cf0, cf1))
#    if (any(cf0 < 0) || any(cf1 < 0)) {
#        warning("Negative prediction: I-value set to NA")
#        I[I < 0 | I > 1] <- NA
#    }
    ## we want model based I (i.e. not mean(Y|z))
    ## thus the 0-1 rescaled version of \beta_1
    ## this is independent from linkinv scaling (and covariates in some sense)
    ## thus comparable across species AND across studies
#    I <- tanh(abs(cf[,2L]))

    out <- data.frame(assoc=h,
        I=tanh(abs(cf[,2L])),
#        I=2*(plogis(abs(cf[,2L]))-0.5),
        null=cfnull,
        mu0=cf0, mu1=cf1,
        beta0=cf[,1L], beta1=cf[,2L],
        logL=ll, logLR=ll-res0$logLik, w=w)
    rownames(out) <- colnames(Z)
    attr(out, "logL_null") <- res0$logLik
#    attr(out, "penalty") <- getOption("ocoptions")$penalty
    attr(out, "H") <- sum(w^2)
    attr(out, "dist") <- if (is.function(dist))
        deparse(substitute(dist)) else dist
    attr(out, "comb") <- Comb
    attr(out, "est") <- Est
    class(out) <- c("opticut1", "data.frame")
    out
}
