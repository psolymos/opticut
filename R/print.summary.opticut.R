print.summary.opticut <- function(x, cut, sort, digits, ...) {
    if (missing(cut))
        cut <- getOption("ocoptions")$cut
    if (missing(digits))
        digits <- max(3L, getOption("digits") - 3L)
    if (missing(sort))
        sort <- getOption("ocoptions")$sort
    sort <- if (is.logical(sort))
        sort[1L] else 1 %in% sort
    xx <- x$summary[, c("split", "assoc", "I", "mu0", "mu1", "logLR", "w")]
    if (sort) {
        xx <- xx[attr(x$bestpart, "row.order"),]
    }
    xx <- xx[xx$logLR >= cut, , drop=FALSE]
    Missing <- nrow(x$summary) - nrow(xx)
    tmp <- if (nrow(xx) > 1L)
        "Best supported models" else "Best supported model"
    TXT <- paste0(tmp, " with logLR >= ",
        format(cut, digits = digits), ":")
    cat("Multivariate opticut results, comb = ", x$comb, ", dist = ", x$dist,
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
