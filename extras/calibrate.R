calibrate <-
function (object, ...)
    UseMethod("calibrate")

calibrate.opticut <-
function(object, ynew, xnew=NULL, cl=NULL, ...)
{
    ynew <- ynew[,colnames(object$Y),drop=FALSE]
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
        #"gaussian:identity"="      y[i,r] ~ dnorm(mu[i,r])",
        "poisson"="      y[i,r] ~ dpois(exp(mu[i,r]))",
        "poisson:log"="      y[i,r] ~ dpois(exp(mu[i,r]))",
        "binomial"="      y[i,r] ~ dbern(ilogit(mu[i,r]))",
        "binomial:logit"="      y[i,r] ~ dbern(ilogit(mu[i,r]))",
        "binomial:cloglog"="      y[i,r] ~ dbern(icloglog(mu[i,r]))",
        "binomial:probit"="      y[i,r] ~ dbern(phi(mu[i,r]))")
    model <- custommodel(unlist(model))
    if (is.null(cl)) {
        jm <- jags.fit(dat, "k", model, ...)
    } else {
        jm <- jags.parfit(cl=cl, dat, "k", model, ...)
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
    list(ynew=ynew,
        xnew=xnew,
        dist=object$dist,
        data=dat,
        model=model,
        niter=nrow(mm),
        gnew=colnames(bp)[gnew],
        pi=t(PI))
}

calibrate.default <-
function(object, ynew, xnew=NULL, cl=NULL, K, ...)
{
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
    model <- custommodel(c(model="model {",
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
        jm <- jags.fit(dat, "k", model, ...)
    } else {
        jm <- jags.parfit(cl=cl, dat, "k", model, ...)
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

