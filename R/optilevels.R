optilevels <-
function(y, x, z = NULL, alpha=0, dist="gaussian", ...)
{
    if (!is.function(dist)) {
        dist <- .opticut_dist(dist, make_dist=TRUE)
        Dist <- strsplit(as.character(dist), ":", fixed=TRUE)[[1L]][1L]
        if (!(Dist %in% c("gaussian", "poisson", "binomial", "negbin",
            "beta", "zip", "zinb")))
            stop("not available for dist=", dist)
    }
    n <- length(y)
    if (is.null(dim(x))) {
        if (length(x) != n)
            stop("length of x must match length of y")
        if (!is.factor(x))
            x <- as.factor(x)
        if (nlevels(x) != length(unique(x)))
            stop("zombie (empty) levels in x")
        X <- model.matrix(~x-1)
        colnames(X) <- levels(x)
    } else {
        if (NROW(x) != n)
            stop("nrows of x must match length of y")
        if (any(colSums(abs(x)) == 0))
            stop("zombie (sum=0) columns in x")
        X <- as.matrix(x)
    }
    if (any(is.na(y)))
        stop("y contains NA")
    if (any(is.na(x)))
        stop("x contains NA")
    if (!is.null(z) && any(is.na(z)))
        stop("z contains NA")
    if (!is.null(z) && NROW(z) != n)
        stop("nrows of z must match length of y")
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
