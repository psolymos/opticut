.summary_opticut <-
function(x, cut=2, sort=TRUE, multi=FALSE)
{
    sort <- if (is.logical(sort))
        sort[1L] else 1 %in% sort
    if (multi) {
        cn <- c("split", "assoc", "I", "G", "J", "logLR")
    } else {
        cn <- c("split", "assoc", "I", "mu0", "mu1", "logLR", "w")
    }
    xx <- x$summary[, cn]
    if (sort) {
        xx <- xx[attr(x$bestpart, "row.order"),]
    }
    xx[xx$logLR >= cut, , drop=FALSE]
}
