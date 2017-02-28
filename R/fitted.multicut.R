fitted.multicut <-
function (object, ...)
{
    if (is.function(object$dist))
        stop("fitted values not available for custom distriutions")
    if (!.opticut_dist(object$dist))
        stop("distribution not recognized")
    Dist <- strsplit(dist, ":", fixed=TRUE)[[1L]][1L]
    Link <- strsplit(dist, ":", fixed=TRUE)[[1L]][2L]
    linkinv <- .opticut1(object$Y[,1L], X=object$X, Z1=NULL,
        dist=object$dist)$linkinv
    Z <- model.matrix(~object$strata)
    K <- ncol(K)
    XX <- cbind(Z, object$X[,-1L,drop=FALSE])
    fun <- function(z, dist, link) {
        cf <- z$coefficients
        if (dist %in% c("zip", "zinb")) {
            if (is.na(link))
                link <- "logit"
            p <- binomial(link)$linkinv(cf[length(cf)])
            return((1 - p) * linkinv(drop(XX %*% cf[-length(cf)])))
        }
        if (dist %in% c("zip2", "zinb2")) {
            p <- linkinv(drop(Z %*% cf[seq_len(K)]))
            lam <- poisson("log")$linkinv(object$X %*% cf[-seq_len(K)])
            return(p * lam)
        }
        if (dist %in% c("zip2", "zinb2")) {
            return(linkinv(drop(XX %*% cf[-length(cf)])))
        }
        linkinv(drop(XX %*% cf))
    }
    fit <- sapply(object$species, f, dist=Dist, link=Link)
    dimnames(fit) <- dimnames(object$Y)
    fit
}
