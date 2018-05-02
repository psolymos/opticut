print.summary.opticut <- function(x, cut, sort, digits, ...) {
    if (missing(cut))
        cut <- getOption("ocoptions")$cut
    if (missing(digits))
        digits <- max(3L, getOption("digits") - 3L)
    if (missing(sort))
        sort <- getOption("ocoptions")$sort

    xx <- .summary_opticut(x, cut=cut, sort=sort)

    Missing <- nrow(x$summary) - nrow(xx)
    tmp <- if (nrow(xx) > 1L)
        "Best supported models" else "Best supported model"
    TXT <- paste0(tmp, " with logLR >= ",
        format(cut, digits = digits), ":")
    cat("Multivariate opticut results, comb = ", x$comb, ", dist = ",
        if (is.function(x$dist)) attr(x$dist, "dist") else x$dist,
        "\n", sep="")
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
            "\n\n", TXT, "\n", sep = "")
    print(format.data.frame(xx, digits=digits), ...)
    cat(x$nsplit, "binary", ifelse(x$nsplit > 1, "splits\n", "split\n"))
    if (Missing)
        cat(Missing, "species not shown\n")
    cat("\n")
    invisible(x)
}
