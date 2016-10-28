## this returns a #habclasses x #bootstrap matrix
bestpart.uncertainty1 <-
function (object, ...)
{
    x <- attr(object, "est")
    if (is.null(x))
        stop("type = 'multi' is required")
    cm <- oComb(x)
    mat <- matrix(0L, nrow(cm), nrow(object))
    rownames(mat) <- rownames(cm)
    for (i in seq_len(nrow(object))) {
        j <- strsplit(as.character(object$best[i]),
            attr(object, "collapse"), fixed=TRUE)[[1L]]
        mat[j,i] <- 1L
    }
#    mat[order(rownames(mat)),]
    mat
}
