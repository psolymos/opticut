print.summary.multicut <- function(x, cut, sort, digits, ...) {
    if (missing(cut))
        cut <- getOption("ocoptions")$cut
    if (missing(digits))
        digits <- max(3L, getOption("digits") - 3L)
    if (missing(sort))
        sort <- getOption("ocoptions")$sort

    xx <- .summary_opticut(x, cut=cut, sort=sort, multi=TRUE)

    Missing <- nrow(x$summary) - nrow(xx)
    tmp <- if (nrow(xx) > 1L)
        "Species models" else "Species model"
    TXT <- paste0(tmp, " with logLR >= ",
        format(cut, digits = digits), ":")
    cat("Multivariate multticut results, dist = ",
        if (is.function(x$dist)) attr(x$dist, "dist") else x$dist,
        "\n", sep="")
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
            "\n\n", TXT, "\n", sep = "")
    print(format(as.data.frame(round(xx, digits)), scientific=FALSE), ...)
    if (Missing)
        cat(Missing, "species not shown\n")
    cat("\n")
    invisible(x)
}
