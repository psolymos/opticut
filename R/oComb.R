## x is a named vector of ranks, referring to factor levels
## in some classification vector, 1=highest abundance.
oComb <-
function(x, collapse)
{
    if (length(x) < 2L)
        stop("length of x must be >1")
    if (missing(collapse))
        collapse <-  getOption("ocoptions")$collapse
    if (is.null(names(x)))
        names(x) <- seq_len(length(x))
    o <- x[order(x, decreasing = FALSE)]
    out <- diag(1L, length(o))
    out[upper.tri(out)] <- 1L
    out <- out[,-ncol(out),drop=FALSE]
    rownames(out) <- names(o)
    colnames(out) <- seq_len(ncol(out))
    for (i in seq_len(ncol(out))) {
        colnames(out)[i] <- paste(names(x)[names(x) %in%
            rownames(out)[out[,i] > 0]],
            collapse = collapse)
    }
    attr(out, "rank") <- o
    out[names(x),]
}
