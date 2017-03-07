.lc_cut <- function(object, fix_fitted=FALSE) {
    n <- table(strata(object))
    mu <- sapply(object$species, "[[", "mu")
    if (fix_fitted)
        mu <- mu + abs(min(mu))
    if (any(mu < 0))
        stop("Negative fitted values found.")
    f <- function(x, n, fix_fitted=FALSE) {
        l <- lorenz(x, n)
        out <- structure(numeric(length(out)), names=names(out))
        out[rownames(l)[-seq_len(which.max(l[,"p"] - l[,"L"]))]] <- 1
        out
    }
    bp <- apply(mu, 2, f, n=n, fix_fitted=fix_fitted)
}
