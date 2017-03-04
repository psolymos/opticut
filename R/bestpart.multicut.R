## this returns an n x #spp matrix
bestpart.multicut <-
function (object, pos_only=FALSE, ...)
{
    fit <- fitted(object)
    lc <- t(apply(fit, 2, function(z) summary(lorenz(z))))
    bp <- t(ifelse(t(fit) > lc[,"x(t)"], 1, 0))
    bp
}
