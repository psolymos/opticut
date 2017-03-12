bsmooth.uncertainty <-
function (object, ...)
{
    if (inherits(object, "uncertainty_opti") && object$type != "multi")
        stop("type='multi' required for bsmooth")
    sapply(lapply(object$uncertainty, bsmooth), rowMeans)
}
