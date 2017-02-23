predict.opticut <-
function (object, gnew=NULL, xnew=NULL, ...)
{
    if (is.null(xnew)) {
        if (!is.null(gnew))
            stop("gnew must be NULL if xnew=NULL")
        fit <- fitted(object)
    } else {
        if (is.null(gnew))
            stop("gnew must be provided")
        if (is.na(object$comb) || object$comb != "rank")
            stop("predict method for new data available only for comb=rank")
        f <- formula(object)
        f[[2]] <- NULL
        bp <- summary(object)$bestpart
        X <- model.matrix(f, xnew)
        linkinv <- .opticut1(object$Y[,1L], X, Z1=NULL, dist=object$dist)$linkinv
        g <- factor(gnew, levels(object$strata))
        if (any(is.na(X)) || any(is.na(g)))
            stop("new data must not have any NA")
        cf <- sapply(names(object$species), function(z) getMLE(object, z)$coef)
        fit <- matrix(0, nrow(X), ncol(cf))
        dimnames(fit) <- list(rownames(X), colnames(cf))
        for (i in seq_len(ncol(cf))) {
            gg <- ifelse(g %in% colnames(bp)[bp[i,] == 1], 1, 0)
            XX <- cbind(X[,1L,drop=FALSE], gg, X[,-1L,drop=FALSE])
            fit[,i] <- linkinv(drop(XX %*% cf[,i]))
        }
    }
    fit
}