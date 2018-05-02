getMLE.opticut <-
function(object, which, vcov=FALSE, ...)
{
    if (missing(which))
        stop("specify which argument")
    if (!length(which))
        stop("which argument must have length 1")
    full_model <- TRUE
    if (!is.function(object$dist)) {
        Dist <- as.character(object$dist)
        Dist <- strsplit(object$dist, ":", fixed=TRUE)[[1L]][1L]
        dist <- object$dist
    } else {
        Dist <- ""
        ## full model cannot be returned for dist=fun with vcov=TRUE
        if (vcov)
            stop("vcov=TRUE cannot be used with custom distribution")
        full_model <- FALSE
        dist <- attr(object$dist, "dist")
    }
    m1 <- .extractOpticut(object, which,
        boot=FALSE,
        internal=TRUE,
        full_model=full_model,
        best=TRUE, ...)[[1L]]
    V <- if (vcov)
        vcov(m1) else NULL
    ## use m1$coef when dist=fun, coef(m1) otherwise
    est <- if (full_model)
        coef(m1) else m1$coef
    # rsf: coef and vcov method returns no intercept
    if (Dist == "rsf") {
        est <- c("(Intercept)"=0, est)
        if (vcov)
            V <- rbind("(Intercept)"=NA, cbind("(Intercept)"=NA, V))
    }
    if (Dist %in% c("zip2", "zinb2")) {
        est <- c(-est[(length(est)-1L):(length(est))],
            est[1:(length(est)-2L)])
        if (vcov)
            V <- V[names(est), names(est)]
    }
    ## these require full model, which is given by full_model
    if (Dist == "negbin")
        attr(est, "theta") <- m1$theta
    if (Dist == "gaussian")
        attr(est, "sigma") <- summary(m1)$sigma
    list(coef=est, vcov=V, dist=dist)
}
