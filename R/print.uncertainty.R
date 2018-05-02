print.uncertainty <-
function(x, ...)
{
    cat("Multivariate opticut uncertainty results",
        ", type = ", x$type, ", B = ", x$B, "\n\n", sep="")
    invisible(x)
}
