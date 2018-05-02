print.summary.uncertainty <-
function(x, sort, digits, ...)
{
    if (missing(digits))
        digits <- max(3L, getOption("digits") - 3L)
    if (missing(sort))
        sort <- getOption("ocoptions")$sort
    inherits(x, "uncertainty_opti")
        cls <- "opticut"
    inherits(x, "uncertainty_multi")
        cls <- "multicut"
    cat("Multivariate ", cls, " uncertainty results",
        "\ntype = ", x$type, ", B = ", x$B,
        ", level = ", format(round(x$level, 2), nsmall=2),
        "\n\n", sep="")
    print(format.data.frame(.summary_uncertainty(x, sort=sort),
        digits=digits), ...)
    invisible(x)
}
