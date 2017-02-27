ipredict.multicut <-
function(object, ynew, xnew=NULL, cl=NULL, K, ...)
{
    if (!all(names(object) == colnames(ynew)))
        stop("names in object and ynew nust match")
    requireNamespace("rjags")
    requireNamespace("dclone")
    ## avoid clashes when running parallel
    ow <- getOption("dcoptions")$overwrite
    dclone::dcoptions("overwrite"=FALSE)
    on.exit(dclone::dcoptions("overwrite"=ow))
    fam <- family(object[[1L]])$family
    link <- family(object[[1L]])$link
    Dist <- paste0(fam, ":", link)
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
    model <- list(
        head=c("model {",
            "  for (i in 1:n) {",
            "    for (r in 1:S) {"),
        dist=NULL,
        mu=NULL,
        mid=c("    }",
            "    k[i] ~ dcat(pi)",
            "  }",
            "  for (r in 1:S) {",
            "    b1[r,1] <- 0",
            "    for (k in 2:K) {",
            "      b1[r,k] <- theta[r,k]",
            "    }",
            "    b0[r] <- theta[r,1]"),
        mvn=NULL,
        tail=c("  }",
            "}"))
    if (!is.null(xnew)) {
        dat$p <- ncol(xnew) - 1L
        dat$x <- xnew[,-1,drop=FALSE] # drop intercept
        model$mu <- "      mu[i,r] <- b0[r] + b1[r,k[i]] + inprod(x[i,], a[r,])"
        model$mvn <- c("    a[r,1:p] <- theta[r,(K+1):(K+p)]",
            "    theta[r,1:(K+p)] ~ dmnorm(cf[1:(K+p),r], prec[1:(K+p),1:(K+p),r])")
    } else {
        model$mu <- "      mu[i,r] <- b0[r] + b1[r,k[i]]"
        model$mvn <- "    theta[r,1:K] ~ dmnorm(cf[1:K,r], prec[1:K,1:K,r])"
    }
    model$dist <- switch(Dist,
        "poisson"="      y[i,r] ~ dpois(exp(mu[i,r]))",
        "poisson:log"="      y[i,r] ~ dpois(exp(mu[i,r]))",
        "binomial"="      y[i,r] ~ dbern(ilogit(mu[i,r]))",
        "binomial:logit"="      y[i,r] ~ dbern(ilogit(mu[i,r]))",
        "binomial:cloglog"="      y[i,r] ~ dbern(icloglog(mu[i,r]))",
        "binomial:probit"="      y[i,r] ~ dbern(phi(mu[i,r]))")
    model <- dclone::custommodel(unlist(model))
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
    Levs <- paste0("Level", seq_len(K))
    rownames(PI) <- Levs
    colnames(PI) <- rownames(ynew)
    gnew <- apply(PI, 2, which.max)
    out <- list(ynew=ynew,
        xnew=xnew,
        dist=Dist,
        data=dat,
        model=as.character(model),
        niter=nrow(mm),
        gnew=factor(Levs[gnew], Levs),
        pi=t(PI))
    class(out) <- c("ipredict.default", "ipredict")
    out
}
