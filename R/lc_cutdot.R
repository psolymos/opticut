.lc_cut <- function(object, fix_fitted=FALSE) {
    n <- table(strata(object))
    mu <- sapply(object$species, "[[", "mu")
    if (fix_fitted)
        mu <- mu + abs(min(mu))
    if (any(mu < 0))
        stop("Negative fitted values found.")
    bp <- apply(mu, 2, .lc_cut1, n=n, fix_fitted=FALSE, bp_only=TRUE)
}
