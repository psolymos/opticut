fitted.opticut <-
function (object, ...)
{
    fit <- sapply(bestmodel(object), fitted)
    rownames(fit) <- rownames(object$Y)
    fit
}
