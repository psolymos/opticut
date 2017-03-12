.lc_cut1 <- function(x, n, fix_fitted=FALSE) {
    if (fix_fitted)
        x <- x + abs(min(x))
    if (any(x < 0))
        stop("Negative fitted values found.")
    l <- lorenz(x, n)
    out <- structure(numeric(length(x)), names=names(x))
    out[rownames(l)[-seq_len(which.max(l[,"p"] - l[,"L"]))]] <- 1
    out
}
