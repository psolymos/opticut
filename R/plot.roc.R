plot.roc <-
function(x, xlab, ylab, ...)
{
    if (missing(xlab))
        xlab <- "False positive rate"
    if (missing(ylab))
        ylab <- "True positive rate"
    plot(x$fpr, x$tpr, type="l", xlab=xlab, ylab=ylab, ...)
    invisible(x)
}
