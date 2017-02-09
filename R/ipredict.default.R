ipredict.default <-
function(object, ynew, xnew=NULL, cl=NULL, K, ...)
{
    requireNamespace("rjags")
    requireNamespace("dclone")
    fam <- family(object[[1L]])$family
    link <- family(object[[1L]])$link
    cf <- sapply(object, coef)
    prec <- array(sapply(object, function(z) solve(vcov(z))),
        dim = c(nrow(cf), nrow(cf), length(object)))
    dat <- list(
        y=ynew,
        pi=rep(1/K, K),
        n=NROW(ynew),
        K=K,
        S=NCOL(ynew),
        cf=cf,
        prec=prec)
    ## this needs to be able to run binomial as well based on family and link
    ## also allow for missing xnew (take ipredict.opticut as example)
    model <- dclone::custommodel(c(model="model {",
        "  for (i in 1:n) {",
        "    for (r in 1:S) {",
        "      y[i,r] ~ dpois(exp(mu[i,r]))", # poisson
        "      mu[i,r] <- b0[r] + b1[r,k[i]] + inprod(x[i,], a[r,])",
        "    }",
        "    k[i] ~ dcat(pi)",
        "  }",
        "  for (r in 1:S) {",
        "    b1[r,1] <- 0",
        "    for (k in 2:K) {",
        "      b1[r,k] <- theta[r,k]",
        "    }",
        "    b0[r] <- theta[r,1]",
        "    a[r,1:p] <- theta[r,(K+1):(K+p)]", # z & x
        "    theta[r,1:(K+p)] ~ dmnorm(cf[1:(K+p),r], prec[1:(K+p),1:(K+p),r])",
        "  }",
        "}"))
    if (!is.null(xnew)) {
        if (is.data.frame(xnew)) {
            ff <- formula(object)
            ff[[2]] <- NULL
            xnew <- model.matrix(ff, xnew)
        }
#        xnew <- xnew[,colnames(object$X),drop=FALSE]
        dat$p <- ncol(xnew) - 1L
        dat$x <- xnew[,-1,drop=FALSE] # drop intercept
    }
    if (is.null(cl)) {
        jm <- dclone::jags.fit(dat, "k", model, ...)
    } else {
        jm <- dclone::jags.parfit(cl=cl, dat, "k", model, ...)
    }
    mm <- as.matrix(jm)
    f <- function(x, K) {
        out <- numeric(K)
        for (k in 1:K)
            out[k] <- sum(x==k)
        out
    }
    PI <- apply(mm, 2, f, K=K) / nrow(mm)
    #rownames(PI) <- colnames(bp)
    colnames(PI) <- rownames(ynew)
    gnew <- apply(PI, 2, which.max)
    list(ynew=ynew,
        xnew=xnew,
#        dist=object$dist,
        data=dat,
        model=model,
        niter=nrow(mm),
        gnew=gnew,
        pi=t(PI))
}
