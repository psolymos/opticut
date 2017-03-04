## this returns an n x #spp matrix
bestpart.multicut <-
function (object, pos_only=FALSE, ...)
{
    fit <- fitted(object)
    lc <- t(apply(fit, 2, function(z) summary(lorenz(z))))
    xt <- lc[,"x(t)"]
    bp <- t(ifelse(t(fit) > xt, 1, 0))
    #attr(bp, "xt") <- xt
    bp
}
