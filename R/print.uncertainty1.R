print.uncertainty1 <-
function(x, ...)
{
    cat("Univariate opticut uncertainty results",
        ", type = ", attr(x, "type"), ", B = ", attr(x, "B"),
        "\n\n", sep="")
    print(summary(x), ...)
    invisible(x)
}
