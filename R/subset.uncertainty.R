subset.uncertainty <-
function(x, subset, ...)
{
    if (any(is.na(subset)))
        stop("subset should not contain NA")
    out <- x
    out$uncertainty <- x$uncertainty[subset]
    if (length(out$uncertainty) < 1L)
        stop("no species left in subset")
    out$Y <- x$Y[,subset,drop=FALSE]
    out
}
