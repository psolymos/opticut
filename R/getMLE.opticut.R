getMLE.opticut <-
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
    } else {
        Dist <- ""
    }
    m1 <- .extractOpticut(object, which,
        boot=FALSE,
        internal=TRUE,
        full_model=vcov,
        #full_model=TRUE,
        best=TRUE, ...)[[1L]]
    if (vcov) {
        est <- coef(m1)
        V <- vcov(m1)
    } else {
        est <- m1$coef
        V <- NULL
    }
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
    if (Dist == "negbin")
        attr(est, "theta") <- m1$theta
    if (Dist == "gaussian")
        attr(est, "sigma") <- summary(m1)$sigma
    list(coef=est, vcov=V, dist=object$dist)
}
