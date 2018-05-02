getMLE.opticut <-
function(object, which, ...)
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
        full_model=TRUE,
        best=TRUE, ...)[[1L]]
    if (Dist == "ordered") {
        est <- c(m1$coefficients, m1$zeta)
        V <- vcov(m1)
        id <- c(length(m1$coefficients)+1,
            1:length(m1$coefficients),
            (length(m1$coefficients)+2):(length(est)))
        est <- est[id]
        V <- V[id, id]
    } else {
        est <- coef(m1)
        V <- vcov(m1)
    }
    if (Dist %in% c("zip2", "zinb2")) {
        est <- est[c((length(est)-1L):(length(est)), 1:(length(est)-2L))]
        V <- V[names(est), names(est)]
    }
    list(coef=est, vcov=V, dist=object$dist)
}
