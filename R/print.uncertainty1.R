print.uncertainty1 <-
function(x, ...)
{
    inherits(x, "uncertainty_opti")
        cls <- "opticut"
    inherits(x, "uncertainty_multi")
        cls <- "multicut"
    cat("Univariate ", cls, " uncertainty results",
        ", type = ", attr(x, "type"), ", B = ", attr(x, "B"),
        "\n\n", sep="")
    print(summary(x), ...)
    invisible(x)
}
