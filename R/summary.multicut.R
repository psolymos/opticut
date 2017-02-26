summary.multicut <- function(object, ...)
{
    object$summary <- t(sapply(object$species, "[[", "mu"))
    object$logLR <- sapply(object$species, "[[", "logLR")
    object$logL_null <- attr(object$species[[1L]], "logL_null")
    object$species <- NULL
    tmp <- t(apply(object$summary, 1, function(z)
        ifelse((z-min(z)) / max(z-min(z)) > 0.5, 1, 0)))
    object$col.order <- order(-colSums(tmp), colnames(tmp))
    object$row.order <- order(ncol(tmp) - rowSums(tmp), decreasing=FALSE)
    class(object) <- c("summary.multicut")
    object
}
