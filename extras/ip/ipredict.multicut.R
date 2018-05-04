ipredict.multicut <-
function(object, ynew, xnew=NULL,
method=c("analytic", "mcmc"), cl=NULL, ...)
{
    if (is.function(object$dist))
        stop("inverse prediction not available for custom distribution")
    dist <- object$dist
    Dist <- .opticut_dist(dist, make_dist=TRUE)
    if (!(Dist %in% c("gaussian", "poisson", "binomial")))
        stop("inverse prediction not available for dist = ", Dist)
    method <- match.arg(method)
    #ynew <- ynew[,colnames(object$Y),drop=FALSE]
    ## new and missing species treated as 0
    sd1 <- setdiff(colnames(ynew), colnames(object$Y))
    sd2 <- setdiff(colnames(object$Y), colnames(ynew))
    if (length(sd1) > 0)
        warning(length(sd1), " in ynew not found in object: dropped")
    if (length(sd2) > 0) {
        stop(length(sd2), " original species not found in ynew: dropped")
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
#    bp <- summary(object)$bestpart
    LEV <- levels(strata(object))
    K <- length(LEV)
    S <- ncol(object$Y)
    mle <- lapply(seq_len(S), function(i)
        getMLE(object, i, vcov=method == "mcmc"))
    cf <- sapply(mle, "[[", "coef")
    if (Dist == "gaussian")
        sig <- sapply(mle, function(z) attr(z$coef, "sigma"))
    if (!is.null(xnew)) {
        if (is.data.frame(xnew)) {
            ff <- formula(object)
            ff[[2]] <- NULL
            xnew <- model.matrix(ff, xnew)
        } else {
            xnew <- data.matrix(xnew)
        }
    }
    if (method == "mcmc") {
        requireNamespace("rjags")
        requireNamespace("dclone")
        ## avoid clashes when running parallel
        ow <- getOption("dcoptions")$overwrite
        dclone::dcoptions("overwrite"=FALSE)
        on.exit(dclone::dcoptions("overwrite"=ow))
        prec <- array(sapply(mle, function(z) solve(z$vcov)),
            dim = c(nrow(cf), nrow(cf), S))
        dat <- list(
            y=ynew,
            pi=rep(1/K, K),
            n=NROW(ynew),
            K=K,
            S=S,
            cf=cf,
            prec=prec)
        if (Dist == "gaussian")
            dat$sigma <- sig
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
        model$dist <- switch(dist,
            "gaussian:identity"="      y[i,r] ~ dnorm(mu[i,r], 1/sigma[r]^2)",
            "gaussian"=         "      y[i,r] ~ dnorm(mu[i,r], 1/sigma[r]^2)",
            "poisson"=          "      y[i,r] ~ dpois(exp(mu[i,r]))",
            "poisson:log"=      "      y[i,r] ~ dpois(exp(mu[i,r]))",
            "binomial"=         "      y[i,r] ~ dbern(ilogit(mu[i,r]))",
            "binomial:logit"=   "      y[i,r] ~ dbern(ilogit(mu[i,r]))",
            "binomial:cloglog"= "      y[i,r] ~ dbern(icloglog(mu[i,r]))",
            "binomial:probit"=  "      y[i,r] ~ dbern(phi(mu[i,r]))")
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
        rownames(PI) <- LEV
        colnames(PI) <- rownames(ynew)
        gnew <- apply(PI, 2, which.max)
        results <- list(
            data=dat,
            model=as.character(model),
            niter=nrow(mm),
            pi=t(PI))
    } else {
        ll <- array(NA, c(dim(ynew), K))
        lls <- matrix(NA, nrow(ynew), K)
        dimnames(ll) <- list(rownames(ynew), colnames(ynew), LEV)
        dimnames(lls) <- list(rownames(ynew), LEV)
        for (i in LEV) {
            pr <- predict(object,
                gnew=factor(rep(i, nrow(ynew)), LEV), xnew=xnew)
            if (Dist == "gaussian") {
                for (j in seq_len(S))
                    ll[,j,i] <- dnorm(ynew[,j], pr[,j], sig[j], log=TRUE)
            }
            if (Dist == "poisson")
                ll[,,i] <- dpois(ynew, pr, log=TRUE)
            if (Dist == "binomial")
                ll[,,i] <- dbinom(ynew, 1, pr, log=TRUE)
            lls[,i] <- rowSums(ll[,,i,drop=FALSE])
        }
        gnew <- apply(lls, 1, which.max)
        results <- list(
            loglik_species=ll,
            loglik=lls)
    }
    out <- list(
        ynew=ynew0,
        xnew=xnew,
        gnew=factor(LEV[gnew], LEV),
        dist=dist,
        method=method,
        results=results)
    class(out) <- c("ipredict.multicut", "ipredict")
    out
}
