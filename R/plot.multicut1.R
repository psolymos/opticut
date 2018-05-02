plot.multicut1 <-
function(x, ylab="Relative abundance", xlab="Strata", ...)
{
    barplot(x$mu, xlab=xlab, ylab=ylab, ...)
    invisible(x)
}

