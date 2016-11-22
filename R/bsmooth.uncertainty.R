bsmooth.uncertainty <-
function (object, ...)
{
    if (object$type != "multi")
        stop("type='multi' required for bsmooth")
    sapply(lapply(object$uncertainty, bsmooth), rowMeans)
}
