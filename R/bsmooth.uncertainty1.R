bsmooth.uncertainty1 <- function (object, ...) {
    if (attr(object,"type") != "multi")
        stop("type='multi' required for bsmooth")
    bp <- bestpart(object)
    t(t(bp) * object$mu1 + (1 - t(bp)) * object$mu0)
}
