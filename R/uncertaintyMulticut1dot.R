## internal wrapper to do single species
.uncertaintyMulticut1 <-
function (object, which=NULL,
type=c("asymp", "boot"), B=99, pb=FALSE, ...)
{
    dots <- setdiff(names(object$call)[-1L],
        c("X", "Y", "formula", "data", "strata", "dist", "comb", "sset", "cl"))
    if (length(dots) > 0)
        stop("Extra arguments detected in opticut call (...)")
    type <- match.arg(type)
    if (missing(which))
        stop("specify which argument")
    if (!length(which))
        stop("which argument must have length 1")
    linkinv <- .opticut1(
        Y=object$Y[,1L],
        X=object$X,
        Z1=NULL,
        dist=object$dist, ...)$linkinv
    scale <- object$scale
    obj <- object$species[[which]]
    #n <- nobs(object)
    K <- length(obj$mu)
    if (type == "asymp") {
        if (length(B) > 1)
            stop("Provide single integer for B.")
        niter <- B
        mle <- getMLE(object, which, vcov=TRUE, ...)
        if (!is.function(object$dist) &&
            .opticut_dist(object$dist, make_dist=TRUE) == "rsf") {
            ## getMLE returns 0 for intercept (NA in vcov)
            cf <- MASS::mvrnorm(niter, mle$coef[-1L],
                mle$vcov[-1L,-1L,drop=FALSE])
            cf <- rbind(mle$coef, cbind(0, cf)) # opticut1
        } else {
            cf <- MASS::mvrnorm(niter, mle$coef, mle$vcov)
            cf <- rbind(mle$coef, cf) # opticut1
        }
    }
    if (type == "boot") {
        if (length(B) == 1) {
            niter <- B
            ## RSF/RSPF requires only used points to be resampled
            if (!is.function(object$dist) &&
                .opticut_dist(object$dist, make_dist=TRUE) %in% c("rsf", "rspf")) {
                avail <- which(object$Y[,1]==0)
                used <- which(object$Y[,1]==1)
                nused <- length(used)
                BB <- replicate(niter, c(sample(used, nused, replace=TRUE), avail))
            } else {
                BB <- replicate(niter, sample.int(n, replace=TRUE))
            }
        } else {
            BB <- B
            niter <- ncol(B)
        }
        nstr <- check_strata(object, BB)
        if (!all(nstr))
            stop("Not all strata represented in resampling")
        m1 <- .extractOpticut(object, which,
            boot=FALSE,
            internal=TRUE,
            full_model=FALSE,
            best=TRUE,
            Z=object$strata, ...)[[1L]]
        cf <- if (pb) {
            t(pbapply::pbapply(BB, 2, function(z, ...) {
                .extractOpticut(object, which,
                    boot=z,
                    internal=TRUE,
                    full_model=FALSE,
                    best=TRUE,
                    Z=object$strata, ...)[[1L]]$coef
            }))
        } else {
            t(apply(BB, 2, function(z, ...) {
                .extractOpticut(object, which,
                    boot=z,
                    internal=TRUE,
                    full_model=FALSE,
                    best=TRUE,
                    Z=object$strata, ...)[[1L]]$coef
            }))
        }
        cf <- rbind(m1$coef, cf)
    }

    mulink <- cf[,seq_len(K),drop=FALSE]
    mulink[,-1] <- mulink[,1] + mulink[,-1,drop=FALSE]
    mu <- linkinv(mulink)
    colnames(mu) <- paste0("mu_", names(obj$mu))
    I <- beta2i(apply(mulink, 1, max) - apply(mulink, 1, min), scale=scale)
    fix <- getOption("ocoptions")$fix_fitted
    if (type == "asymp")
        bp <- apply(mu, 1, function(z)
            .lc_cut1(x=structure(z, names=names(obj$mu)),
                n=table(strata(object)), fix_fitted=fix))
    if (type == "boot")
        bp <- cbind(.lc_cut1(x=structure(mu[1L,], names=names(obj$mu)),
                n=table(object$strata), fix_fitted=fix),
            sapply(seq_len(niter), function(i)
                .lc_cut1(x=structure(mu[i+1L,], names=names(obj$mu)),
                n=table(object$strata[BB[,i]]),fix_fitted=fix)))
    lab1 <- character(ncol(bp))
    for (i in seq_len(ncol(bp))) {
        lab1[i] <- paste(rownames(bp)[bp[,i] == 1],
            collapse=getOption("ocoptions")$collapse)
    }
    out <- data.frame(best=lab1, I=I, mu)
    class(out) <- c("uncertainty1_multi", "uncertainty1", "data.frame")
    attr(out, "B") <- niter
    attr(out, "type") <- type
    attr(out, "scale") <- scale
    attr(out, "collapse") <- object$collapse
    out
}
