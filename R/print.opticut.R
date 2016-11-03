print.opticut <- function(x, digits, ...) {
    if (missing(digits))
        digits <- max(3L, getOption("digits") - 3L)
    cat("Multivariate opticut results, comb = ", x$comb, ", dist = ",
        if (is.function(x$dist)) attr(x$dist, "dist") else x$dist,
        "\n", sep="")
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
            "\n\n", sep = "")
    cat(length(x$species), "species, ")
    cat(x$nsplit, ifelse(x$nsplit > 1, "binary splits\n", "binary split\n"))
    cat("\n")
    invisible(x)
}

