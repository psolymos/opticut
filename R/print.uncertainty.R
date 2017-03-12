print.uncertainty <-
function(x, ...)
{
    inherits(object, "uncertainty_opti")
        cls <- "opticut"
    inherits(object, "uncertainty_multi")
        cls <- "multicut"
    cat("Multivariate", cls, "uncertainty results",
        ", type = ", x$type, ", B = ", x$B, "\n\n", sep="")
    invisible(x)
}
