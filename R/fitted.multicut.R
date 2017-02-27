fitted.multicut <-
function (object, ...)
{
    linkinv <- .opticut1(object$Y[,1L], X=object$X, Z1=NULL,
        dist=object$dist)$linkinv
    fit <- sapply(object$species, function(z) {
        linkinv(drop(object$X %*% z$coefficients))
    })
    dimnames(fit) <- dimnames(object$Y)
    fit
}
