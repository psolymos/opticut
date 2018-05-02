.lc_cut <- function(x, n, fix_fitted=FALSE) {
    if (fix_fitted && any(x < 0))
        x <- x + abs(min(x))
    if (any(x < 0))
        stop("Negative fitted values found.")
    l <- lorenz(x, n)
    bp <- structure(numeric(length(x)), names=names(x))
    bp[rownames(l)[-seq_len(which.max(l[,"p"] - l[,"L"]))]] <- 1
    bp
}
