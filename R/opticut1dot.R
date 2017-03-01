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
        dist <- strsplit(dist, ":", fixed=TRUE)[[1L]]
        if (length(dist) > 1L) {
            link <- dist[2L]
            dist <- dist[1L]
        } else {
            link <- NULL
        }
        dist <- match.arg(dist, .opticut_dist())
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
            if (dist == "gaussian")
                attr(cf, "sigma") <- summary(mod)$sigma
            ll <- as.numeric(logLik(mod))
            linv <- family(mod)$linkinv
        }
        if (dist == "negbin") {
            mod <- MASS::glm.nb(Y ~ .-1, data=XX, ...)
            cf <- coef(mod)
            attr(cf, "theta") <- mod$theta
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
            ## log link used for count based contrasts
            linv <- poisson("log")$linkinv
        }
## zip2 & zinb2 implementation:
## - MLE returns unmodified coefs (P of 0 in ZI)
## - .opticut1 returns:
##       -1*coef[1:2]
##       linkinv: binomial(link)$linkinv(eta)
## - asymp uncertainty uses MLE, thus have to invert and use linkinv after
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
            ## only symmetric function can be simply inverted on link scale
            if (!(link %in% c("logit", "probit")))
                stop("only logit and probit link allowed for zip2 and zinb2")
            mod <- pscl::zeroinfl(Y ~ X-1 | ZZ-1, dist=Dist,
                link=link, ...)
            ## logit/probit/etc used for (1-phi) based contrasts
            linv <- binomial(link)$linkinv
            cf <- c(-coef(mod, "zero"), coef(mod, "count"))
            ll <- as.numeric(logLik(mod))
        }
        if (dist %in% c("rsf", "rspf")) {
            ## formula interface used (rsf, rspf, and not rsf.fit)
            ## thus sanity checks are made in ResourceSelection
            ## regarding covariates/identifiability.
            ## Note: m=0 is provided, otherwise resampling is difficult,
            ## leave m != 0 cases to customization (dist=fun)
            if (dist == "rsf") {
                if (!is.null(link) && link != "log")
                    warning("link argument ignored for dist='rsf' (log link used by default)")
                link <- "log" # this is needed for linkinv below
                ## --- rsf implementation uses PL of Lele 2009 ---
                ## note: the call uses link='log', no need to provide link=link here!
                ## intercept is not reported by rsf
                ## and this can cause problems in X %*% theta
                mod <- if (ncol(XX) < 2L) {
                    ResourceSelection::rsf.null(Y, m=0, ...)
                } else {
                    ResourceSelection::rsf(Y ~ .,
                        data=XX[,-1,drop=FALSE], m=0, B=0, ...)
                }
                cf <- mod$results$par
                cf[1L] <- 0
                ## --- glm implementation uses ML ---
                ## can be used for null model as well
#                mod <- stats::glm(Y ~ .-1, data=XX,
#                    family=binomial("logit"), ...)
#                cf <- coef(mod)
#                ll <- as.numeric(logLik(mod))
                ## link is still log !!!
            } else {
                if (is.null(link))
                    link <- "logit" # this is needed for linkinv below
                mod <- ResourceSelection::rspf(Y ~ ., data=XX[,-1,drop=FALSE],
                    link=link, m=0, B=0, ...)
                cf <- mod$coefficients
#                ll <- as.numeric(mod$loglik)
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
