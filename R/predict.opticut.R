predict.opticut <-
function (object, gnew=NULL, xnew=NULL, ...)
{
    if (is.null(xnew) && is.null(gnew)) {
        fit <- fitted(object)
    } else {
        if (is.na(object$comb) || object$comb != "rank")
            stop("predict method for new data available only for comb=rank")
        if (is.null(gnew))
            stop("gnew must be provided")
        if (is.null(xnew)) {
            xnew <- matrix(1, length(gnew), 1L)
        }
        bp <- summary(object)$bestpart
        ff <- formula(object)
        ff[[2]] <- NULL
        X <- if (is.data.frame(xnew))
            model.matrix(ff, xnew) else data.matrix(xnew)
        Dist <- strsplit(object$dist, ":", fixed=TRUE)[[1L]][1L]
        Link <- strsplit(object$dist, ":", fixed=TRUE)[[1L]][2L]
        linkinv <- .opticut1(object$Y[,1L], X=object$X, Z1=NULL,
            dist=object$dist)$linkinv
        g <- factor(gnew, levels(object$strata))
        if (any(is.na(X)) || any(is.na(g)))
            stop("new data must not have any NA")
        if (nrow(X) != length(g))
            stop("length of gnew must equal nrow(xnew)")
        cf <- sapply(names(object$species), function(z) getMLE(object, z)$coef)
        fit <- matrix(0, nrow(X), ncol(cf))
        dimnames(fit) <- list(rownames(X), colnames(cf))
        for (i in seq_len(ncol(cf))) {
            gg <- ifelse(g %in% colnames(bp)[bp[i,] == 1], 1, 0)
            fit[,i] <- .predict_dist(cf[,i],
                dist=Dist, link=Link, X=X, Z=gg, linkinv=linkinv)
#            XX <- cbind(X[,1L,drop=FALSE], gg, X[,-1L,drop=FALSE])
#            fit[,i] <- linkinv(drop(XX %*% cf[,i]))
        }
    }
    fit
}
