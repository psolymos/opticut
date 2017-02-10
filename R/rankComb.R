rankComb <-
function(Y, X, Z, dist="gaussian", collapse, ...)
{
    if (missing(X))
        X <- matrix(1L, length(Y), 1L)
    if (!is.factor(Z))
        stop("Z must be a factor")
    if (length(Y) != length(Z))
        stop("length(Y) must equal length(Z)")
    if (missing(collapse))
        collapse <-  getOption("ocoptions")$collapse
    Z0 <- model.matrix(~Z, model.frame(~Z, na.action=na.pass))
    m <- .opticut1(Y, X, Z1=Z0[,-1L,drop=FALSE],
        linkinv=TRUE, dist=dist, ...)
    lc <- c(m$coef[1], m$coef[1] + m$coef[2:ncol(Z0)])
    names(lc) <- levels(Z)
    x <- rank(-lc)
    oc <- oComb(x, collapse = collapse)
    out <- oc[match(Z, rownames(oc)),,drop=FALSE]
    attr(out, "est") <- m$linkinv(lc)
    attr(out, "collapse") <- collapse
    attr(out, "comb") <- "rank"
    out
}
