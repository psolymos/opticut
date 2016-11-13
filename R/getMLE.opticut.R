getMLE.opticut <-
function(object, which, ...)
{
    if (missing(which))
        stop("specify which argument")
    if (!length(which))
        stop("which argument must have length 1")
        ## full model cannot be returned for dist=fun
    m1 <- .extractOpticut(object, which,
        boot=FALSE,
        internal=TRUE,
        full_model=TRUE,
        best=TRUE, ...)[[1L]]
    if (object$dist == "ordered") {
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
    list(coef=est, vcov=V, dist=object$dist)
}
