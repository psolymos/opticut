fitted.opticut <-
function (object, ...)
{
    if (is.function(object$dist))
        stop("fitted values not available for custom distriutions")
    if (!.opticut_dist(object$dist))
        stop("distribution not recognized")
    fit <- sapply(bestmodel(object), fitted)
    rownames(fit) <- rownames(object$Y)
    fit
}
