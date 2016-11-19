print.summary.lorenz <-
function(x, digits, ...)
{
    if (missing(digits))
        digits <- max(3L, getOption("digits") - 3L)
    cat("Lorenz curve summary\n\n")
    print(x[-1], digits = digits)
    cat("\n")
    invisible(x)
}
