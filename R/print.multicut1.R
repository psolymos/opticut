print.multicut1 <- function(x, digits, ...)
{
    if (missing(digits))
        digits <- max(3L, getOption("digits") - 3L)
    cat("Univariate multicut results, dist = ", attr(x, "dist"),
        "\nlogLR = ", format(x$logLR, digits = digits),
        " (logL_null = ", format(attr(x, "logL_null"), digits = digits),
        ")\n\nExpected values:\n", sep="")
    print(x$mu, digits=digits, ...)
    cat("\n")
    invisible(x)
}
