print.uncertainty1 <-
function(x, ...)
{
    inherits(x, "uncertainty1_opti")
        cls <- "opticut"
    inherits(x, "uncertainty1_multi")
        cls <- "multicut"
    cat("Univariate ", cls, " uncertainty results",
        ", type = ", attr(x, "type"), ", B = ", attr(x, "B"),
        "\n\n", sep="")
    print(summary(x), ...)
    invisible(x)
}
