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
            ## formula interface used (rsf, rspf, and not rspf.fit)
            ## thus sanity checks are made in ResourceSelection
            ## regarding covariates/identifiability.
            ## Note: m=0 is provided, otherwise resampling is difficult,
            ## leave m not 0 cases to customization (dist=fun)
            if (dist == "rsf") {
                if (!is.null(link) && link != "log")
                    warning("link argument ignored for dist='rsf' (log link used by default)")
                link <- "log" # this is needed for linkinv below
                ## note: the call uses link='log', no need to provide link=link here!
                mod <- ResourceSelection::rsf(Y ~ .,
                    data=XX[,-1,drop=FALSE], m=0, B=0, ...)
                ## intercept is not reported by rsf
                ## and this can cause problems in X %*% theta
                cf <- c(0, mod$coefficients)
            } else {
                if (is.null(link))
                    link <- "logit" # this is needed for linkinv below
                mod <- ResourceSelection::rspf(Y ~ ., data=XX[,-1,drop=FALSE],
                    link=link, m=0, B=0, ...)
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
