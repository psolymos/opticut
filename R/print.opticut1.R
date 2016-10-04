print.opticut1 <- function(x, cut, sort, digits, ...)
{
    if (missing(cut))
        cut <- getOption("ocoptions")$cut
    if (missing(digits))
        digits <- max(3L, getOption("digits") - 3L)
    if (missing(sort))
        sort <- getOption("ocoptions")$sort
    ## sorting only rows/species here
    sort <- if (is.logical(sort))
        sort[1L] else 1 %in% sort
    xx <- x
    xx$assoc <- .parseAssoc(xx)
    xx <- xx[, c("assoc", "I", "mu0", "mu1", "logLR", "w")]
    if (sort)
        xx <- xx[order(xx$logLR, decreasing=TRUE),]
    if (any(xx$logLR >= cut)) {
        SHOW <- which(xx$logLR >= cut)
        tmp <- if (length(SHOW) > 1L)
            "Best supported models" else "Best supported model"
        TXT <- paste0(tmp, " with logLR >= ",
            format(cut, digits = digits), ":")
    } else {
        SHOW <- 1L
        TXT <- paste0("Best supported model:")
    }
    xx <- xx[SHOW,,drop=FALSE]
    cat("Univariate opticut results, comb = ", attr(x, "comb"),
        ", dist = ", attr(x, "dist"),
        "\nI = ", format(xx[1L,"I"], digits = digits),
        "; w = ", format(xx[1L,"w"], digits = digits),
        "; H = ", format(attr(x, "H"), digits = digits),
        "; logL_null = ", format(attr(x, "logL_null"), digits = digits),
        "\n\n", TXT, "\n", sep="")
    print.data.frame(xx, digits=digits, ...)
    DROP <- nrow(x) - nrow(xx)
    if (DROP > 0) {
        cat(nrow(x), " binary ",
            ifelse(nrow(x) > 1, "splits", "split"),
            " (", DROP,
            ifelse(DROP > 1, " models", " model"),
            " not shown)\n", sep="")
    } else {
        cat(nrow(x), "binary",
            ifelse(nrow(x) > 1, "splits", "split"), "\n")
    }
    cat("\n")
    invisible(x)
}
