bsmooth.uncertainty1 <- function (object, ...) {
    if (inherits(object, "uncertainty_opti")) {
        if (attr(object,"type") != "multi")
            stop("type='multi' required for bsmooth")
        bp <- bestpart(object)
        out <- t(t(bp) * object$mu1 + (1 - t(bp)) * object$mu0)
    }
    if (inherits(object, "uncertainty_opti")) {
        stop("bsmooth is not available for multicut based objects.")
    }
    out
}
