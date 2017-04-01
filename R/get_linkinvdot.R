.get_linkinv <- function(object, ...) {
    .opticut1(object$Y[,1L], dist=object$dist, ...)$linkinv
}
