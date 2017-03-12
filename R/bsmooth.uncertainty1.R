bsmooth.uncertainty1 <- function (object, ...) {
    if (inherits(object, "uncertainty1_opti"))
        stop("bsmooth is not available for multicut based objects.")
    if (inherits(object, "uncertainty1_opti")) {
        if (attr(object,"type") != "multi")
            stop("type='multi' required for bsmooth")
        bp <- bestpart(object)
        return(t(t(bp) * object$mu1 + (1 - t(bp)) * object$mu0))
    }
}
