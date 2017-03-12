print.uncertainty1 <-
function(x, ...)
{
    inherits(object, "uncertainty_opti")
        cls <- "opticut"
    inherits(object, "uncertainty_multi")
        cls <- "multicut"
    cat("Univariate", cls, "uncertainty results",
        ", type = ", attr(x, "type"), ", B = ", attr(x, "B"),
        "\n\n", sep="")
    print(summary(x), ...)
    invisible(x)
}
