.lc_cut1 <- function(x, n, fix_fitted=FALSE, bp_only=TRUE) {
    if (fix_fitted)
        x <- x + abs(min(x))
    if (any(x < 0))
        stop("Negative fitted values found.")
    l <- lorenz(x, n)
    bp <- structure(numeric(length(x)), names=names(x))
    bp[rownames(l)[-seq_len(which.max(l[,"p"] - l[,"L"]))]] <- 1
    if (bp_only) {
        out <- bp
    } else {
        attr(l, "bp") <- bp
        mu0 <- sum((1-bp)*x*n) / sum((1-bp)*n)
        mu1 <- sum(bp*x*n) / sum(bp*n)
        mu <- c(mu0=mu0, mu1=mu1)
        attr(l, "mu") <- mu
        out <- l
    }
    out
}
