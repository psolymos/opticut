## do not cut, but sort???
print.summary.uncertainty <-
function(x, digits, ...)
{
    if (missing(digits))
        digits <- max(3L, getOption("digits") - 3L)
    cat("Multivariate opticut uncertainty results",
        "\ntype = ", x$type, ", B = ", x$B,
        ", level = ", format(round(x$level, 2), nsmall=2),
        "\n\n", sep="")
    uct <- x$uctab
    uct <- uct[order(uct$split, -uct$R, -uct$I),]
    print(format.data.frame(uct, digits=digits), ...)
    invisible(x)
}
