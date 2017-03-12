print.uncertainty <-
function(x, ...)
{
    inherits(x, "uncertainty_opti")
        cls <- "opticut"
    inherits(x, "uncertainty_multi")
        cls <- "multicut"
    cat("Multivariate ", cls, " uncertainty results",
        ", type = ", x$type, ", B = ", x$B, "\n\n", sep="")
    invisible(x)
}
