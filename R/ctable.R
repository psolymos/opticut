ctable <-
function(x, y)
{
    if (length(y) != length(x))
        stop("Hey! x and y lengths must match")
    x <- as.factor(x)
    x <- droplevels(x)
    levs <- sort(levels(x))
    x <- factor(x, levels=levs)
    y <- factor(y, levels=levs)
    if (any(is.na(y)))
        stop("NAs, and levels in y that are not in x are not allowed, sorry")
    as.matrix(table(Reference=x, Prediction=y))
}
