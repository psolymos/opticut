subset.multicut <-
function(x, subset, ...)
{
    if (any(is.na(subset)))
        stop("subset should not contain NA")
    out <- x
    out$species <- x$species[subset]
    if (length(out$species) < 1L)
        stop("no species left in subset")
    out$Y <- x$Y[,subset,drop=FALSE]
    out
}
