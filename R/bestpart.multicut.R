## this returns an n x #spp matrix
bestpart.multicut <-
function (object, ...)
{
    bp <- sapply(object$species, "[[", "bestpart")
    bp[match(strata(object), rownames(bp)),,drop=FALSE]
}
