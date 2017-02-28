fitted.multicut <-
function (object, ...)
{
    if (is.function(object$dist))
        stop("fitted values not available for custom distriutions")
    if (!.opticut_dist(object$dist))
        stop("distribution not recognized")
    Dist <- strsplit(object$dist, ":", fixed=TRUE)[[1L]][1L]
    Link <- strsplit(object$dist, ":", fixed=TRUE)[[1L]][2L]
    linkinv <- .opticut1(object$Y[,1L], X=object$X, Z1=NULL,
        dist=object$dist)$linkinv
    Z <- model.matrix(~object$strata)[,-1L,drop=FALSE]
    fit <- sapply(object$species, function(z)
        .predict_dist(z$coefficient,
            dist=Dist, link=Link, X=object$X, Z=Z, linkinv=linkinv))
    dimnames(fit) <- dimnames(object$Y)
    fit
}
