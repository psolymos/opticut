bsmooth.uncertainty1 <- function (object, ...) {
    bp <- bestpart(object)
    t(t(bp) * object$mu1 + (1 - t(bp)) * object$mu0)
}
