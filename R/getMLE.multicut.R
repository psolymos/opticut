getMLE.multicut <-
function(object, which, vcov=FALSE, ...)
{
    if (missing(which))
        stop("specify which argument")
    if (!length(which))
        stop("which argument must have length 1")
    ## full model cannot be returned for dist=fun
    if (!is.function(object$dist)) {
        Dist <- as.character(object$dist)
        Dist <- strsplit(object$dist, ":", fixed=TRUE)[[1L]][1L]
        dist <- object$dist
    } else {
        Dist <- ""
        dist <- attr(object$dist, "dist")
    }
    K <- nlevels(object$strata)
    if (vcov)
        m1 <- bestmodel(object, which, ...)[[1L]]
    est <- if (vcov)
        coef(m1) else object$species[[which]]$coefficients
    V <- if (vcov)
        vcov(m1) else NULL
    if (vcov && Dist == "rsf") {
        est <- c("(Intercept)"=0, est)
        V <- rbind("(Intercept)"=NA, cbind("(Intercept)"=NA, V))
    }
    if (vcov && Dist %in% c("zip2", "zinb2")) {
        est <- est[c((length(est)-(K-1L)):(length(est)), 1:(length(est)-K))]
        V <- V[names(est), names(est)]
    }
    if (Dist == "negbin")
        attr(est, "theta") <- if (vcov)
            attr(est, "theta") else m1$theta
    if (Dist == "gaussian")
        attr(est, "sigma") <- if (vcov)
            summary(m1)$sigma else attr(est, "sigma")
    list(coef=est, vcov=V, dist=dist)
}
