optilevels <-
function(y, x, z = NULL, alpha=0, dist="gaussian", ...)
{
    if (is.null(dim(x))) {
        if (!is.factor(x))
            x <- as.factor(x)
        if (nlevels(x) != length(unique(x)))
            stop("zombie (empty) levels in x")
        na.op <- options(na.action="na.pass")
        on.exit(options(na.action=na.op))
        X <- model.matrix(~x-1)
        options(na.action=na.op)
        colnames(X) <- levels(x)
    } else {
        if (any(colSums(abs(x)) == 0))
            stop("zombie (sum=0) columns in x")
        X <- as.matrix(x)
    }
    out <- .optilevels(Y=y, X=X, Z=z, alpha=alpha, dist=dist, ...)
    levs <- list()
    for (i in seq_len(length(out$ranklist))) {
        levi <- sapply(1:max(out$rank[i,]), function(j)
            paste(colnames(out$coef)[out$rank[i,] == j],
            collapse=getOption("ocoptions")$collapse))
        levs[[i]] <- levi[out$rank[i,]]
        names(levs[[i]]) <- colnames(out$coef)
    }
    out$levels <- levs
    out$factor <- is.factor(x)
    out$call <- match.call()
    class(out) <- "optilevels"
    out
}
