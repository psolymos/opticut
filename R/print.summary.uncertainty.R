## do not cut, but sort???
print.summary.uncertainty <-
function(x, ...)
{
    cat("Multivariate opticut uncertainty results",
        ", type = ", x$type, ", B = ", x$B,
        "\nlevel = ", format(round(x$level, 2), nsmall=2),
        "\n\n", sep="")
    uct <- x$uctab
    #uct <- uct[order(uct$split, -uct$PC),]
    print(uct, ...)
    invisible(x)
}
