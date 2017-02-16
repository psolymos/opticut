ipredict.opticut <-
function(object, ynew, xnew=NULL, cl=NULL, ...)
{
    requireNamespace("rjags")
    requireNamespace("dclone")
    ## avoid clashes when running parallel
    ow <- getOption("dcoptions")$overwrite
    dclone::dcoptions("overwrite"=FALSE)
    on.exit(dclone::dcoptions("overwrite"=ow))
    #ynew <- ynew[,colnames(object$Y),drop=FALSE]
    ## new and missing species treated as 0
    sd1 <- setdiff(colnames(ynew), colnames(object$Y))
    sd2 <- setdiff(colnames(object$Y), colnames(ynew))
    if (length(s1) > 0)
        warning(length(s1), " in ynew not found in object: dropped")
    if (length(s2) > 0) {
        stop(length(s2), " original species not found in ynew: dropped")
        object <- subset(object,
            colnames(object$Y)[colnames(object$Y) %in% colnames(ynew)])
    }
    ynew0 <- ynew
    ynew <- matrix(0, nrow(ynew), ncol(object$Y))
    dimnames(ynew) <- list(rownames(ynew0), colnames(object$Y))
    cn <- intersect(colnames(ynew0), colnames(object$Y))
    ynew[,cn] <- ynew0[,cn]
    if (!is.null(xnew) && ncol(object$X) < 2L)
        xnew <- NULL
    bp <- summary(object)$bestpart
    K <- ncol(bp)
    S <- nrow(bp)
    mle <- lapply(seq_len(S), function(i) getMLE(object, i))
    cf <- sapply(mle, "[[", "coef")
    prec <- array(sapply(mle, function(z) solve(z$vcov)),
        dim = c(nrow(cf), nrow(cf), S))
    dat <- list(
        y=ynew,
        bp=bp,
        pi=rep(1/K, K),
        n=NROW(ynew),
        S=S,
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
            "    b0[r] <- theta[r,1]",
            "    b1[r] <- theta[r,2]"),
        mvn=NULL,
        tail=c(    "  }",
            "}"))
    if (!is.null(xnew)) {
        if (is.data.frame(xnew)) {
            ff <- formula(object)
            ff[[2]] <- NULL
            xnew <- model.matrix(ff, xnew)
        }
        xnew <- xnew[,colnames(object$X),drop=FALSE]
        dat$p <- ncol(object$X) - 1L
        dat$x <- xnew[,-1,drop=FALSE] # drop intercept
        model$mu <- "      mu[i,r] <- b0[r] + bp[r,k[i]]*b1[r] + inprod(x[i,], a[r,])"
        model$mvn <- c("    a[r,1:p] <- theta[r,3:(2+p)]",
            "    theta[r,1:(2+p)] ~ dmnorm(cf[1:(2+p),r], prec[1:(2+p),1:(2+p),r])")
    } else {
        model$mu <- "      mu[i,r] <- b0[r] + bp[r,k[i]]*b1[r]"
        model$mvn <- "    theta[r,1:2] ~ dmnorm(cf[1:2,r], prec[1:2,1:2,r])"
    }
    model$dist <- switch(object$dist,
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
    rownames(PI) <- colnames(bp)
    colnames(PI) <- rownames(ynew)
    gnew <- apply(PI, 2, which.max)
    out <- list(ynew=ynew0,
        xnew=xnew,
        dist=object$dist,
        data=dat,
        model=as.character(model),
        niter=nrow(mm),
        gnew=factor(colnames(bp)[gnew], colnames(bp)),
        pi=t(PI))
    class(out) <- c("ipredict.opticut", "ipredict")
    out
}
