## Z1 is:
## * NULL (this is used to fit model for H0, i.e. no partition)
## * a single column from Z matrix
## * or a matrix itself (design matrix w/o intercept)
.opticut1 <-
function(Y, X, Z1=NULL,
dist="gaussian", linkinv, full_model=FALSE, ...)
{
    if (missing(linkinv))
        linkinv <- is.null(Z1)
    if (is.null(Z1)) {
        XX <- as.data.frame(X)
    } else {
        if (is.null(dim(X)))
            X <- data.matrix(X)
        XX <- as.data.frame(cbind(X[,1,drop=FALSE], Z1, X[,-1,drop=FALSE]))
    }
    if (!is.function(dist)) {
        dist <- as.character(dist)
        dist <- strsplit(dist, ":", fixed=TRUE)[[1]]
        if (length(dist) > 1L) {
            link <- dist[2L]
            dist <- dist[1L]
        } else {
            link <- NULL
        }
        dist <- match.arg(dist,
            c("gaussian","poisson","binomial","negbin",
            "beta","zip","zinb","ordered", "rsf", "rspf",
            "zip2", "zinb2"))
        if (dist %in% c("gaussian", "poisson", "binomial")) {
            if (is.null(link))
                link <- switch(dist,
                    "gaussian"="identity",
                    "poisson"="log",
                    "binomial"="logit")
            Family <- switch(dist,
                "gaussian"=gaussian(link),
                "poisson"=poisson(link),
                "binomial"=binomial(link))
            mod <- stats::glm(Y ~ .-1, data=XX, family=Family, ...)
            cf <- coef(mod)
            ll <- as.numeric(logLik(mod))
            linv <- family(mod)$linkinv
        }
        if (dist == "negbin") {
            mod <- MASS::glm.nb(Y ~ .-1, data=XX, ...)
            cf <- coef(mod)
            ll <- as.numeric(logLik(mod))
            linv <- family(mod)$linkinv
        }
        if (dist == "beta") {
            mod <- betareg::betareg(Y ~ .-1, data=XX, ...)
            cf <- coef(mod)
            ll <- as.numeric(logLik(mod))
            linv <- mod$link$mean$linkinv
        }
        if (dist %in% c("zip", "zinb")) {
            if (is.null(link))
                link <- "logit"
            Dist <- switch(dist,
                "zip"="poisson",
                "zinb"="negbin")
            mod <- pscl::zeroinfl(Y ~ .-1 | 1, data=XX, dist=Dist,
                link=link, ...)
            cf <- coef(mod)
            ll <- as.numeric(logLik(mod))
            linv <- mod$linkinv
        }
        if (dist %in% c("zip2", "zinb2")) {
            ZZ <- if (is.null(Z1)) {
                data.matrix(X[,1,drop=FALSE])
            } else {
                data.matrix(cbind(X[,1,drop=FALSE], Z1))
            }
            if (is.null(link))
                link <- "logit"
            Dist <- switch(dist,
                "zip2"="poisson",
                "zinb2"="negbin")
            mod <- pscl::zeroinfl(Y ~ X-1 | ZZ-1, dist=Dist,
                link=link, ...)
            cf <- c(coef(mod, "zero"), coef(mod, "count"))
            ll <- as.numeric(logLik(mod))
            linv <- function(eta) 1 - binomial(link)$linkinv(eta)
        }
        if (dist == "ordered") {
            if (!is.null(list(...)$method))
                if (list(...)$method != "logistic")
                    stop("Sorry but only logisic model allowed for ordered.")
            if (!is.null(link))
                if (link != "logistic")
                    stop("Sorry but only logisic model allowed for ordered.")
            Y <- as.ordered(Y)
            if (nlevels(Y) > 2) { # ordinal
                ## need to keep the intercept
                if (ncol(XX)==1) {
                    if (!is.null(list(...)$data))
                        stop("Note: data argument should not be provided as part of ...")
                    mod <- MASS::polr(Y ~ 1, method="logistic", ...)
                } else {
                    mod <- MASS::polr(Y ~ ., data=data.frame(XX[,-1,drop=FALSE]),
                        method="logistic", ...)
                }
                cf <- c(mod$zeta[1], coef(mod))
            } else {
                mod <- stats::glm(Y ~ .-1, data=XX, family=binomial("logit"), ...)
                cf <- coef(mod)
            }
            ll <- as.numeric(logLik(mod))
            linv <- binomial("logit")$linkinv
        }
        if (dist %in% c("rsf", "rspf")) {
            m <- list(...)$m
            if (dist == "rsf") {
                if (is.null(m))
                    stop("'m' must be provided, have you checked help('rsf') ?")
                if (!is.null(link) && link != "log")
                    warning("link argument ignored for dist='rsf' (log link used by default)")
                link <- "log" # this is needed for linkinv below
                ## note: the call uses link='log', no need to provide link=link here!
                mod <- ResourceSelection::rsf(Y ~ ., data=XX[,-1,drop=FALSE], ...)
                ## intercept is not reported by rsf
                ## and this can cause problems in X %*% theta
                cf <- c(0, mod$coefficients)
            } else {
                if (is.null(m))
                    stop("'m' must be provided, have you checked help('rspf') ?")
                if (is.null(link))
                    link <- "logit" # this is needed for linkinv below
                mod <- ResourceSelection::rspf(Y ~ ., data=XX[,-1,drop=FALSE],
                    link=link, ...)
                cf <- mod$coefficients
            }
            ll <- as.numeric(mod$loglik)
            linv <- binomial(link)$linkinv
        }
        if (!linkinv)
            linv <- NULL
        out <- list(coef=cf, logLik=ll, linkinv=linv)
    } else {
        if (full_model)
            stop("Unable to return full model for custom distribution function.")
        out <- dist(Y, XX, linkinv, ...)
    }
    if (full_model)
        mod else out
}

## todo: if Z inherits from class optilevels,
## use that as best binary partition (check no. of levels)

## Y is abundance vector
## X is model matrix for nuisance variables
## Z is design matrix for binary splits or a factor (using rankComb)
opticut1 <-
function(Y, X, Z, dist="gaussian", sset=NULL, ...)
{
    if (missing(X))
        X <- matrix(1, length(Y), 1)
    X <- data.matrix(X)
    if (is.null(rownames(X)))
        rownames(X) <- seq_len(nrow(X))
    if (is.factor(Z)) {
        Z <- rankComb(Y, X, Z, dist=dist, ...)
        Est <- attr(Z, "est")
        Comb <- "rank"
    } else {
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
        stop("Dimension mismatch: check you input.")
    if (is.null(rownames(Z))) {
        warning("Row names added to binary split matrix Z (it was NULL). You are welcome.")
        rownames(Z) <- apply(Z, 1, paste, collapse="")
    }
    N <- ncol(Z)
    if (!is.null(sset)) {
        Y <- Y[sset]
        X <- X[sset,,drop=FALSE]
        Z <- Z[sset,,drop=FALSE]
    }
    res0 <- .opticut1(Y, X, Z1=NULL, dist=dist, ...)
    cf <- matrix(0, N, length(res0$coef)+1)
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

## this is the main user interface
opticut <-
function(formula, data, strata, dist="gaussian",
comb=c("rank", "all"), sset=NULL, cl=NULL, ...)
{
    comb <- match.arg(comb)
    if (missing(data))
        data <- parent.frame()
    if (missing(strata))
        stop("It looks like that strata is missing.")
    Strata <- deparse(substitute(strata))
    if (Strata %in% names(data))
        strata <- data[[Strata]]
    mf <- match.call(expand.dots = FALSE)
    mm <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, mm)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    Y <- model.response(mf, "numeric")
    Y <- data.matrix(Y)
    if (is.null(colnames(Y)))
        colnames(Y) <- paste("Species", seq_len(ncol(Y)))
    if (any(duplicated(colnames(Y))))
        warning("Duplicate column names found and renamed in LHS.")
    colnames(Y) <- make.names(colnames(Y))
    ff <- formula
    ff[[2]] <- NULL
    mt <- terms(ff, data = data)
    X <- model.matrix(mt, mf)

    if (is.null(dim(strata))) {
        if (nchar(getOption("ocoptions")$collapse) < 1)
            stop("nchar(getOption('ocoptions')$collapse) must be > 0")
        if (is.ordered(strata)) {
            warning("ordering in strata ignored")
            class(strata) <- "factor"
        }
        strata <- droplevels(as.factor(strata)) # factor
        ## make syntactically valid names
        #levels(strata) <- make.names(levels(strata), unique = TRUE)
        ## make sure that collapse is not in levels
        if (any(grepl(getOption('ocoptions')$collapse, levels(strata), fixed=TRUE)))
            stop("Collapse value found in levels.")
        if (comb == "rank") {
            Z <- strata
        }
        if (comb == "all") {
            Z <- allComb(strata) # matrix
        }
    } else {
        Z <- as.matrix(strata) # matrix
        if (getOption("ocoptions")$check_comb && !checkComb(Z))
            stop("Guess what! Complementary design variables found:\nuse 'checkComb'")
        #colnames(Z) <- make.names(colnames(Z), unique = TRUE)
        comb <- NA # user supplied matrix, not checked
    }
    Y <- data.matrix(Y)

    if (!is.function(dist)) {
        Dist <- strsplit(as.character(dist), ":", fixed=TRUE)[[1]][1]
        Dist <- match.arg(Dist,
            c("gaussian","poisson","binomial","negbin",
            "beta","zip","zinb","ordered", "rsf", "rspf",
            "zip2", "zinb2"))
        ## sanity check for rsf/rspf
        if (Dist %in% c("rsf", "rspf")) {
            if (ncol(Y) > 1L)
                stop("rsf/rspf is only available for single species in RHS")
            if (identical(as.character(ff[[2]]), "1"))
                stop("invalid formula, no covariates")
            factonly <- all(unique(sapply(mf, .MFclass)[-1]) %in% c("ordered", "factor"))
            if (factonly && dist == "rspf")
                stop("provide at least 1 continuous covariate for RSPF")
        }
        ## sanity check for rsf/rspf
        if (Dist == "ordered") {
            if (ncol(Y) > 1L)
                stop("ordered is only available for single species in RHS")
        }
    }

    ## sequential
    if (is.null(cl)) {
        ## show progress bar
        if (!ncol(Y) < 2L) {
            pbo <- pboptions(type = "none")
            on.exit(pboptions(pbo))
        }
        if (getOption("ocoptions")$try_error) {
            res <- pbapply::pbapply(Y, 2, function(yy, ...)
                try(opticut1(Y=yy, X=X, Z=Z, dist=dist, sset=sset, ...)), ...)
        } else {
            res <- pbapply::pbapply(Y, 2, function(yy, ...)
                opticut1(Y=yy, X=X, Z=Z, dist=dist, sset=sset, ...), ...)
        }
    ## parallel
    } else {
        ## snow type cluster
        if (inherits(cl, "cluster")) {
            if (length(cl) < 2)
                stop("Are you kidding? Set cl to utilize at least 2 workers.")
            parallel::clusterEvalQ(cl, library(opticut))
            .oc_envir <- new.env()
            assign("dist", dist, envir=.oc_envir)
            assign("X", X, envir=.oc_envir)
            assign("Z", X, envir=.oc_envir)
            assign("Y", Y, envir=.oc_envir)
            assign("sset", sset, envir=.oc_envir)
#            parallel::clusterExport(cl, c("Y", "X","Z","dist"), envir=e)
            parallel::clusterExport(cl, ".oc_envir")
            if (getOption("ocoptions")$try_error) {
                res <- parallel::parLapply(cl, seq_len(ncol(Y)), function(i, ...)
                    try(opticut1(Y=.oc_envir$Y[,i],
                        X=.oc_envir$X, Z=.oc_envir$Z,
                        dist=.oc_envir$dist, sset=.oc_envir$sset, ...)), ...)
            } else {
                res <- parallel::parLapply(cl, seq_len(ncol(Y)), function(i, ...)
                    opticut1(Y=.oc_envir$Y[,i],
                        X=.oc_envir$X, Z=.oc_envir$Z,
                        dist=.oc_envir$dist, sset=.oc_envir$sset, ...), ...)
            }
            #parallel::clusterEvalQ(cl, rm(list=c("Y", "X","Z","dist")))
            parallel::clusterEvalQ(cl, rm(list=".oc_envir"))
            parallel::clusterEvalQ(cl, detach(package:opticut))
        ## forking
        } else {
            if (.Platform$OS.type == "windows" && cl != 1)
                stop("Unfortunately forking (cl > 1) does not work on Windows:",
                     "try cl as a cluster instead.")
            if (cl < 2)
                stop("Are you kidding? Set cl to utilize at least 2 workers.")
            if (getOption("ocoptions")$try_error) {
                res <- parallel::mclapply(seq_len(ncol(Y)), function(i, ...)
                    try(opticut1(Y=Y[,i], X=X, Z=Z, dist=dist, sset=sset, ...)),
                    mc.cores = cl, ...)
            } else {
                res <- parallel::mclapply(seq_len(ncol(Y)), function(i, ...)
                    opticut1(Y=Y[,i], X=X, Z=Z, dist=dist, sset=sset, ...),
                    mc.cores = cl, ...)
            }
        }
    }
    if (getOption("ocoptions")$try_error) {
        Failed <- sapply(res, inherits, "try-error")
        failed <- names(res)[Failed]
        if (any(Failed)) {
            if (length(failed) == length(res))
                stop("Bad news: opticut failed for all species.")
            warning("Bad news: opticut failed for ", length(failed),
                " out of ", length(res), " species.")
        }
    } else {
        Failed <- logical(length(res))
        failed <- character(0)
    }
    NOBS <- if (is.null(sset))
        NROW(Y) else NROW(data.matrix(Y)[sset,,drop=FALSE])
    out <- list(call=match.call(),
        species=res[!Failed],
        X=X,
        Y=Y[,!Failed,drop=FALSE],
        strata=Z,
        nobs=NOBS,
        sset=sset,
        nsplit=if (is.factor(Z)) # strata as factor implies K-1 splits
            (nlevels(Z) - 1L) else ncol(Z),
        collapse=getOption("ocoptions")$collapse,
        dist=if (is.function(dist))
            deparse(substitute(dist)) else dist,
        comb=comb,
        failed=failed,
        collapse=getOption("ocoptions")$collapse)
    class(out) <- "opticut"
    out
}
