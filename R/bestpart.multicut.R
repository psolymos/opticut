## this returns an n x #spp matrix
bestpart.multicut <-
function (object, ...)
{
    fit <- fitted(object)
    if (getOption("ocoptions")$fix_fitted)
        fit <- fit + abs(min(fit))
    if (any(fit < 0))
        stop("Negative fitted values found.")
    lc <- t(apply(fit, 2, function(z) summary(lorenz(z))))
    out <- t(ifelse(t(fit) >= lc[,"x(t)"], 1, 0))
    #attr(out, "p") <- lc[,"p(t)"]
    rownames(out) <- strata(object)
    out
}
