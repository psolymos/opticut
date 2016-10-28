.summary_opticut <- function(x, cut=2, sort=TRUE) {

    sort <- if (is.logical(sort))
        sort[1L] else 1 %in% sort
    xx <- x$summary[, c("split", "assoc", "I", "mu0", "mu1", "logLR", "w")]
    if (sort) {
        xx <- xx[attr(x$bestpart, "row.order"),]
    }
    xx[xx$logLR >= cut, , drop=FALSE]
}
