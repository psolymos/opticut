.get_linkinv <- function(object, ...) {
    .opticut1(object$Y[,1L], X=object$X, Z1=NULL, dist=object$dist, ...)$linkinv
}
