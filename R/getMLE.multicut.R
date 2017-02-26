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
    } else {
        Dist <- ""
    }
    K <- nlevels(object$strata)
    if (vcov)
        m1 <- .extractOpticut(object, which,
            boot=FALSE,
            internal=TRUE,
            full_model=TRUE,
            best=TRUE,
            Z=object$strata, ...)[[1L]]
    if (Dist == "ordered") {
        est <- c(m1$coefficients, m1$zeta)
        V <- vcov(m1)
        id <- c(length(m1$coefficients)+1,
            1:length(m1$coefficients),
            (length(m1$coefficients)+2):(length(est)))
        est <- est[id]
        V <- V[id, id]
    } else {
        est <- if (vcov)
            coef(m1) else object$species[[which]]$coefficients
        V <- if (vcov)
            vcov(m1) else NULL
    }
    if (Dist %in% c("zip2", "zinb2")) {
        est <- est[c((length(est)-(K-1L)):(length(est)), 1:(length(est)-K))]
        if (vcov)
            V <- V[names(est), names(est)]
    }
    list(coef=est, vcov=V, dist=object$dist)
}
