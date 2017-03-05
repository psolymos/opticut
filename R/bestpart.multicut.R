## this returns an n x #spp matrix
bestpart.multicut <-
function (object, ...)
{
    fit <- fitted(object)
    lc <- t(apply(fit, 2, function(z) summary(lorenz(z))))
    xt <- lc[,"x(t)"]
    out <- t(ifelse(t(fit) > xt, 1, 0))
    #attr(bp, "xt") <- xt
    rownames(out) <- strata(object)
    out
}
